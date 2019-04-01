/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

// Overload of uBLAS prod function with MKL/GSL implementations
#include "votca/xtp/qmpackage.h"
#include "votca/xtp/aomatrix.h"
#include <boost/algorithm/string.hpp>

namespace votca {
namespace xtp {
using std::flush;
void QMPackage::ReorderOutput(Orbitals& orbitals) {
  BasisSet dftbasisset;
  dftbasisset.LoadBasisSet(_basisset_name);
  if (!orbitals.hasQMAtoms()) {
    throw std::runtime_error("Orbitals object has no QMAtoms");
  }

  AOBasis dftbasis;
  dftbasis.AOBasisFill(dftbasisset, orbitals.QMAtoms());
  // necessary to update nuclear charges on qmatoms
  if (_write_pseudopotentials) {
    BasisSet ecps;
    ecps.LoadPseudopotentialSet(_ecp_name);
    AOBasis ecpbasis;
    ecpbasis.ECPFill(ecps, orbitals.QMAtoms());
  }

  if (orbitals.hasAOOverlap()) {
    dftbasis.ReorderMatrix(orbitals.AOOverlap(), getPackageName(), "xtp");
    CTP_LOG(ctp::logDEBUG, *_pLog) << "Reordered Overlap matrix" << flush;
  }
  if (orbitals.hasAOVxc()) {
    dftbasis.ReorderMatrix(orbitals.AOVxc(), getPackageName(), "xtp");
    CTP_LOG(ctp::logDEBUG, *_pLog) << "Reordered VXC matrix" << flush;
  }
  if (orbitals.hasMOCoefficients()) {
    dftbasis.ReorderMOs(orbitals.MOCoefficients(), getPackageName(), "xtp");
    CTP_LOG(ctp::logDEBUG, *_pLog) << "Reordered MOs" << flush;
  }

  return;
}

void QMPackage::ReorderMOsBack(Orbitals& orbitals) {
  BasisSet dftbasisset;
  dftbasisset.LoadBasisSet(_basisset_name);
  if (!orbitals.hasQMAtoms()) {
    throw std::runtime_error("Orbitals object has no QMAtoms");
  }
  AOBasis dftbasis;
  dftbasis.AOBasisFill(dftbasisset, orbitals.QMAtoms());
  dftbasis.ReorderMOs(orbitals.MOCoefficients(), "xtp", getPackageName());
  return;
}

std::vector<std::vector<double> > QMPackage::SplitMultipoles(
    ctp::APolarSite* aps) {

  std::vector<std::vector<double> > multipoles_split;

  const tools::vec pos = aps->getPos() * tools::conv::nm2ang;
  tools::vec tot_dpl = tools::vec(0.0);
  if (_with_polarization) {
    tot_dpl += aps->getU1();
  }
  if (aps->getRank() > 0) {
    tot_dpl += aps->getQ1();
  }
  // Calculate virtual charge positions
  double a = _dpl_spacing;      // this is in nm
  double mag_d = abs(tot_dpl);  // this is in e * nm
  if (mag_d > 1e-9) {
    tools::vec dir_d = tot_dpl.normalize();
    tools::vec A =
        pos + 0.5 * a * dir_d * tools::conv::nm2ang;  // converted to AA
    tools::vec B = pos - 0.5 * a * dir_d * tools::conv::nm2ang;
    double qA = mag_d / a;
    double qB = -qA;
    multipoles_split.push_back({A.getX(), A.getY(), A.getZ(), qA});
    multipoles_split.push_back({B.getX(), B.getY(), B.getZ(), qB});
  }

  if (aps->getRank() > 1) {
    tools::matrix components = aps->getQ2cartesian();
    tools::matrix::eigensystem_t system;
    components.SolveEigensystem(system);
    double a = 2 * _dpl_spacing;
    for (unsigned i = 0; i < 3; i++) {

      double q = system.eigenvalues[i] / (a * a);
      if (std::abs(q) < 1e-9) {
        continue;
      }
      tools::vec vec1 =
          pos + 0.5 * a * system.eigenvecs[i] * tools::conv::nm2ang;
      tools::vec vec2 =
          pos - 0.5 * a * system.eigenvecs[i] * tools::conv::nm2ang;

      multipoles_split.push_back({vec1.getX(), vec1.getY(), vec1.getZ(), q});
      multipoles_split.push_back({vec2.getX(), vec2.getY(), vec2.getZ(), q});
    }
  }
  return multipoles_split;
}

void QMPackage::setMultipoleBackground(
    std::vector<std::shared_ptr<ctp::PolarSeg> > PolarSegments) {
  if (PolarSegments.size() == 0) {
    std::cout << "WARNING::The Multipole Background has no entries!"
              << std::endl;
    return;
  }
  _PolarSegments = PolarSegments;
  _write_charges = true;

  WriteChargeOption();
}

std::vector<std::string> QMPackage::GetLineAndSplit(
    std::ifstream& input_file, const std::string separators) {
  std::string line;
  getline(input_file, line);
  boost::trim(line);
  std::vector<std::string> row;
  boost::algorithm::split(row, line, boost::is_any_of(separators),
                          boost::algorithm::token_compress_on);
  return row;
}

std::vector<std::string> QMPackage::FindUniqueElements(
    const std::vector<QMAtom*> atoms) {
  std::vector<std::string> result;
  for (QMAtom* atom : atoms) {
    bool exists = false;
    if (result.size() == 0) {
      exists = false;
    } else {
      for (const std::string& type : result) {
        if (atom->getType() == type) {
          exists = true;
          break;
        }
      }
    }
    if (!exists) {
      result.push_back(atom->getType());
    }
  }
  return result;
}

}  // namespace xtp
}  // namespace votca
