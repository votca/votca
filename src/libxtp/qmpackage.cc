/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <boost/algorithm/string.hpp>
#include <votca/xtp/ecpaobasis.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmpackage.h>
#include <votca/xtp/qmpackagefactory.h>

namespace votca {
namespace xtp {
using std::flush;

tools::Property QMPackage::ParseCommonOptions(const tools::Property& options) {

  std::string key = "package";
  // std::string key = "";

  _settings.read_property(options, key);
  char* votca_share = getenv("VOTCASHARE");
  if (votca_share == nullptr) {
    std::cout << "Warning: VOTCASHARE environment variable not defined\n";
  } else {
    Settings qmpackage_defaults{key};
    qmpackage_defaults.load_from_xml(this->FindDefaultsFile());
    _settings.amend(qmpackage_defaults);
  }
  _settings.validate();

  _charge = _settings.get<Index>("charge");
  _spin = _settings.get<Index>("spin");
  _basisset_name = _settings.get("basisset");

  if (_settings.has_key("cleanup")) {
    _cleanup = _settings.get("cleanup");
  }

  if (getPackageName() != "xtp") {
    _scratch_dir = _settings.get("scratch");
  }
  return _settings.to_property("package");
}

void QMPackage::ReorderOutput(Orbitals& orbitals) const {
  if (!orbitals.hasQMAtoms()) {
    throw std::runtime_error("Orbitals object has no QMAtoms");
  }

  AOBasis dftbasis = orbitals.SetupDftBasis();
  // necessary to update nuclear charges on qmatoms
  if (orbitals.hasECPName()) {
    ECPBasisSet ecps;
    ecps.Load(orbitals.getECPName());
    ECPAOBasis ecpbasis;
    ecpbasis.Fill(ecps, orbitals.QMAtoms());
  }

  if (orbitals.hasMOs()) {
    ReorderMOsToXTP(orbitals.MOs().eigenvectors(), dftbasis);
    XTP_LOG(Log::info, *_pLog) << "Reordered MOs" << flush;
  }

  return;
}

Eigen::MatrixXd QMPackage::ReorderMOsBack(const Orbitals& orbitals) const {
  if (!orbitals.hasQMAtoms()) {
    throw std::runtime_error("Orbitals object has no QMAtoms");
  }
  AOBasis dftbasis = orbitals.SetupDftBasis();
  Eigen::MatrixXd result = orbitals.MOs().eigenvectors();
  ReorderMOsToNative(result, dftbasis);
  return result;
}

std::vector<QMPackage::MinimalMMCharge> QMPackage::SplitMultipoles(
    const StaticSite& aps) const {

  std::vector<QMPackage::MinimalMMCharge> multipoles_split;
  // Calculate virtual charge positions
  double a = _settings.get<double>("dipole_spacing");  // this is in a0
  double mag_d = aps.getDipole().norm();               // this is in e * a0
  const Eigen::Vector3d dir_d = aps.getDipole().normalized();
  const Eigen::Vector3d A = aps.getPos() + 0.5 * a * dir_d;
  const Eigen::Vector3d B = aps.getPos() - 0.5 * a * dir_d;
  double qA = mag_d / a;
  double qB = -qA;
  if (std::abs(qA) > 1e-12) {
    multipoles_split.push_back(MinimalMMCharge(A, qA));
    multipoles_split.push_back(MinimalMMCharge(B, qB));
  }

  if (aps.getRank() > 1) {
    const Eigen::Matrix3d components = aps.CalculateCartesianMultipole();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.computeDirect(components);
    double a2 = 2 * a;
    for (Index i = 0; i < 3; i++) {
      double q = es.eigenvalues()[i] / (a2 * a2);
      const Eigen::Vector3d vec1 =
          aps.getPos() + 0.5 * a2 * es.eigenvectors().col(i);
      const Eigen::Vector3d vec2 =
          aps.getPos() - 0.5 * a2 * es.eigenvectors().col(i);
      multipoles_split.push_back(MinimalMMCharge(vec1, q));
      multipoles_split.push_back(MinimalMMCharge(vec2, q));
    }
  }
  return multipoles_split;
}

std::vector<std::string> QMPackage::GetLineAndSplit(
    std::ifstream& input_file, const std::string separators) const {
  std::string line;
  getline(input_file, line);
  boost::trim(line);
  tools::Tokenizer tok(line, separators);
  return tok.ToVector();
}

std::string QMPackage::FindDefaultsFile() const {
  auto xmlFile = std::string(getenv("VOTCASHARE")) +
                 std::string("/xtp/xml/qmpackage_defaults.xml");

  return xmlFile;
}

void QMPackage::ReorderMOs(Eigen::MatrixXd& v,
                           const std::vector<Index>& order) const {
  // Sanity check
  if (v.rows() != Index(order.size())) {
    throw std::runtime_error("Size mismatch in ReorderMOs " +
                             std::to_string(v.rows()) + ":" +
                             std::to_string(order.size()));
  }
  // actual swapping of coefficients
  for (Index s = 1, d; s < (Index)order.size(); ++s) {
    for (d = order[s]; d < s; d = order[d]) {
      ;
    }
    if (d == s) {
      while (d = order[d], d != s) {
        v.row(s).swap(v.row(d));
      }
    }
  }
}

void QMPackage::ReorderMOsToXTP(Eigen::MatrixXd& v,
                                const AOBasis& basis) const {

  std::vector<Index> multiplier = getMultiplierVector(basis);

  std::vector<Index> order = getReorderVector(basis);
  ReorderMOs(v, order);
  MultiplyMOs(v, multiplier);
  return;
}

void QMPackage::ReorderMOsToNative(Eigen::MatrixXd& v,
                                   const AOBasis& basis) const {

  std::vector<Index> multiplier = getMultiplierVector(basis);
  MultiplyMOs(v, multiplier);
  std::vector<Index> order = getReorderVector(basis);
  std::vector<Index> reverseorder = invertOrder(order);
  ReorderMOs(v, reverseorder);
  return;
}

void QMPackage::MultiplyMOs(Eigen::MatrixXd& v,
                            const std::vector<Index>& multiplier) const {
  // Sanity check
  if (v.cols() != Index(multiplier.size())) {
    std::cerr << "Size mismatch in MultiplyMOs" << v.cols() << ":"
              << multiplier.size() << std::endl;
    throw std::runtime_error("Abort!");
  }
  for (Index i = 0; i < v.cols(); i++) {
    v.row(i) = multiplier[i] * v.row(i);
  }
  return;
}

std::vector<Index> QMPackage::getMultiplierVector(const AOBasis& basis) const {
  std::vector<Index> multiplier;
  multiplier.reserve(basis.AOBasisSize());
  // go through basisset
  for (const AOShell& shell : basis) {
    std::vector<Index> shellmultiplier = getMultiplierShell(shell);
    multiplier.insert(multiplier.end(), shellmultiplier.begin(),
                      shellmultiplier.end());
  }
  return multiplier;
}

std::vector<Index> QMPackage::getMultiplierShell(const AOShell& shell) const {
  // multipliers were all found using code, hard to establish
  std::vector<Index> multiplier;
  for (Index i = 0; i < shell.getNumFunc(); i++) {
    multiplier.push_back(ShellMulitplier()[shell.getOffset() + i]);
  }
  return multiplier;
}

std::vector<Index> QMPackage::getReorderVector(const AOBasis& basis) const {
  std::vector<Index> reorder;
  reorder.reserve(basis.AOBasisSize());
  // go through basisset
  for (const AOShell& shell : basis) {
    std::vector<Index> shellreorder = getReorderShell(shell);
    std::for_each(shellreorder.begin(), shellreorder.end(),
                  [&reorder](Index& i) { i += Index(reorder.size()); });
    reorder.insert(reorder.end(), shellreorder.begin(), shellreorder.end());
  }
  return reorder;
}

std::vector<Index> QMPackage::getReorderShell(const AOShell& shell) const {
  std::vector<Index> reorder;
  for (Index i = 0; i < shell.getNumFunc(); i++) {
    // Shellreorder tells how much a certain funktion has to be shifted with
    // reference to votca ordering
    reorder.push_back(ShellReorder()[shell.getOffset() + i] + i);
  }
  return reorder;
}

std::vector<Index> QMPackage::invertOrder(
    const std::vector<Index>& order) const {

  std::vector<Index> neworder = std::vector<Index>(order.size());
  for (Index i = 0; i < Index(order.size()); i++) {
    neworder[order[i]] = Index(i);
  }
  return neworder;
}

}  // namespace xtp
}  // namespace votca
