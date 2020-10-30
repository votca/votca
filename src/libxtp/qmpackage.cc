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

// Third party includes
#include <boost/algorithm/string.hpp>

// Local VOTCA includes
#include "votca/tools/globals.h"
#include "votca/xtp/ecpaobasis.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmpackage.h"
#include "votca/xtp/qmpackagefactory.h"

namespace votca {
namespace xtp {
using std::flush;

tools::Property QMPackage::ParseCommonOptions(const tools::Property& options) {

  std::string key = "package";

  _settings.read_property(options, key);

  if (tools::VotcaShareSet()) {
    Settings qmpackage_defaults{key};
    qmpackage_defaults.load_from_xml(this->FindDefaultsFile());
    _settings.amend(qmpackage_defaults);
  } else {
    std::cout << "Warning: VOTCASHARE environment variable not defined\n";
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
    OrbReorder reorder(ShellReorder(), ShellMulitplier());
    reorder.reorderOrbitals(orbitals.MOs().eigenvectors(), dftbasis);
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
  bool reverseOrder = true;
  OrbReorder reorder(ShellReorder(), ShellMulitplier(), reverseOrder);
  reorder.reorderOrbitals(result, dftbasis);
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
  return tools::GetVotcaShare() + "/xtp/data/qmpackage_defaults.xml";
}

}  // namespace xtp
}  // namespace votca
