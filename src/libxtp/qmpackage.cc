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

namespace votca {
namespace xtp {
using std::flush;

void QMPackage::ParseCommonOptions(tools::Property& options) {

  std::string key = "package";

  _settings.read_property(options, key);
  Settings qmpackage_template{key};
  qmpackage_template.load_from_xml(this->FindTemplateFile());
  _settings.merge(qmpackage_template);
  _settings.validate();

  _charge = _settings.get<Index>("charge");
  _spin = _settings.get<Index>("spin");
  _basisset_name = _settings.get("basisset");
  _auxbasisset_name = _settings.get("auxbasisset");

  _cleanup = _settings.get("cleanup");
  _write_guess = _settings.get<bool>("read_guess");
  _write_charges = _settings.get<bool>("write_charges");

  if (getPackageName() != "xtp") {
    _executable = _settings.get("executable");
    _scratch_dir = _settings.get("scratch");
    // TODO: REMOVE THIS SECTION AFTER CLEANING ALL THE PACKAGES
    if (_settings.exists("options")) {
      _settings.get("option");
    }
  }
  if (_settings.exists("ecp")) {
    _write_pseudopotentials = true;
    _ecp_name = _settings.get("ecp");
  }
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
    dftbasis.ReorderMOs(orbitals.MOs().eigenvectors(), getPackageName(), "xtp");
    XTP_LOG(Log::info, *_pLog) << "Reordered MOs" << flush;
  }

  return;
}

Eigen::MatrixXd QMPackage::ReorderMOsBack(const Orbitals& orbitals) const {
  BasisSet dftbasisset;
  dftbasisset.Load(_basisset_name);
  if (!orbitals.hasQMAtoms()) {
    throw std::runtime_error("Orbitals object has no QMAtoms");
  }
  AOBasis dftbasis;
  dftbasis.Fill(dftbasisset, orbitals.QMAtoms());
  Eigen::MatrixXd result = orbitals.MOs().eigenvectors();
  dftbasis.ReorderMOs(result, "xtp", getPackageName());
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
    double a2 = 2 * _settings.get<double>("dipole_spacing");
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

std::string QMPackage::FindTemplateFile() const {
  auto xmlFile = std::string(getenv("VOTCASHARE")) +
                 std::string("/xtp/packages/qmpackage_template.xml");

  return xmlFile;
}

}  // namespace xtp
}  // namespace votca
