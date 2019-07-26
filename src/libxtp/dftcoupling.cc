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

#include <votca/xtp/aomatrix.h>
#include <votca/xtp/dftcoupling.h>

#include <boost/format.hpp>
#include <boost/progress.hpp>
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {

using boost::format;
using std::flush;

void DFTcoupling::Initialize(tools::Property& options) {

  std::string key = "";
  _degeneracy =
      options.ifExistsReturnElseReturnDefault<bool>(key + "degeneracy", 0.0);
  _numberofstatesA =
      options.ifExistsReturnElseReturnDefault<int>(key + "levA", 1);
  _numberofstatesB =
      options.ifExistsReturnElseReturnDefault<int>(key + "levB", 1);
}

void DFTcoupling::WriteToProperty(tools::Property& type_summary,
                                  const Orbitals& orbitalsA,
                                  const Orbitals& orbitalsB, int a, int b) {
  double J = getCouplingElement(a, b, orbitalsA, orbitalsB);
  tools::Property& coupling = type_summary.add("coupling", "");
  double energyA = orbitalsA.getMOEnergy(a) * tools::conv::hrt2ev;
  double energyB = orbitalsB.getMOEnergy(b) * tools::conv::hrt2ev;
  coupling.setAttribute("levelA", a);
  coupling.setAttribute("levelB", b);
  coupling.setAttribute("j", (format("%1$1.6e") % J).str());
  coupling.setAttribute("eA", (format("%1$1.6e") % energyA).str());
  coupling.setAttribute("eB", (format("%1$1.6e") % energyB).str());
}

void DFTcoupling::Addoutput(tools::Property& type_summary,
                            const Orbitals& orbitalsA,
                            const Orbitals& orbitalsB) {
  tools::Property& dftcoupling = type_summary.add(Identify(), "");
  dftcoupling.setAttribute("homoA", orbitalsA.getHomo());
  dftcoupling.setAttribute("homoB", orbitalsB.getHomo());
  tools::Property& hole_summary = dftcoupling.add("hole", "");
  // hole hole
  for (int a = Range_orbA.first; a <= orbitalsA.getHomo(); ++a) {
    for (int b = Range_orbB.first; b <= orbitalsB.getHomo(); ++b) {
      WriteToProperty(hole_summary, orbitalsA, orbitalsB, a, b);
    }
  }
  tools::Property& electron_summary = dftcoupling.add("electron", "");
  // electron-//electron
  for (int a = orbitalsA.getLumo(); a <= Range_orbA.first + Range_orbA.second;
       ++a) {
    for (int b = orbitalsB.getLumo(); b <= Range_orbB.first + Range_orbB.second;
         ++b) {
      WriteToProperty(electron_summary, orbitalsA, orbitalsB, a, b);
    }
  }
  return;
}

std::pair<int, int> DFTcoupling::DetermineRangeOfStates(
    const Orbitals& orbital, int numberofstates) const {
  const Eigen::VectorXd& MOEnergies = orbital.MOs().eigenvalues();
  if (std::abs(MOEnergies(orbital.getHomo()) - MOEnergies(orbital.getLumo())) <
      _degeneracy) {
    throw std::runtime_error(
        "Homo Lumo Gap is smaller than degeneracy. "
        "Either your degeneracy is too large or your Homo and Lumo are "
        "degenerate");
  }

  int minimal = orbital.getHomo() - numberofstates + 1;
  int maximal = orbital.getLumo() + numberofstates - 1;

  std::vector<int> deg_min = orbital.CheckDegeneracy(minimal, _degeneracy);
  for (int i : deg_min) {
    if (i < minimal) {
      minimal = i;
    }
  }

  std::vector<int> deg_max = orbital.CheckDegeneracy(maximal, _degeneracy);
  for (int i : deg_max) {
    if (i > maximal) {
      maximal = i;
    }
  }

  std::pair<int, int> result;
  result.first = minimal;                 // start
  result.second = maximal - minimal + 1;  // size

  return result;
}

double DFTcoupling::getCouplingElement(int levelA, int levelB,
                                       const Orbitals& orbitalsA,
                                       const Orbitals& orbitalsB) const {

  int levelsA = Range_orbA.second;

  if (_degeneracy != 0) {
    std::vector<int> list_levelsA =
        orbitalsA.CheckDegeneracy(levelA, _degeneracy);
    std::vector<int> list_levelsB =
        orbitalsA.CheckDegeneracy(levelB, _degeneracy);

    double JAB_sq = 0;

    for (int iA : list_levelsA) {
      for (int iB : list_levelsB) {
        double JAB_one_level = JAB(iA - 1, iB - 1 + levelsA);
        JAB_sq += JAB_one_level * JAB_one_level;
      }
    }
    return std::sqrt(JAB_sq / (list_levelsA.size() * list_levelsB.size())) *
           tools::conv::hrt2ev;
  } else {
    return JAB(levelA - 1, levelB - 1 + levelsA) * tools::conv::hrt2ev;
  }
}

/**
 * \brief evaluates electronic couplings
 * @param _orbitalsA molecular orbitals of molecule A
 * @param _orbitalsB molecular orbitals of molecule B
 * @param _orbitalsAB molecular orbitals of the dimer AB
 */
void DFTcoupling::CalculateCouplings(const Orbitals& orbitalsA,
                                     const Orbitals& orbitalsB,
                                     const Orbitals& orbitalsAB) {

  XTP_LOG(logDEBUG, *_pLog) << "Calculating electronic couplings" << flush;

  CheckAtomCoordinates(orbitalsA, orbitalsB, orbitalsAB);

  // constructing the direct product orbA x orbB
  int basisA = orbitalsA.getBasisSetSize();
  int basisB = orbitalsB.getBasisSetSize();

  if ((basisA == 0) || (basisB == 0)) {
    throw std::runtime_error("Basis set size is not stored in monomers");
  }

  Range_orbA = DetermineRangeOfStates(orbitalsA, _numberofstatesA);
  Range_orbB = DetermineRangeOfStates(orbitalsB, _numberofstatesB);

  int levelsA = Range_orbA.second;
  int levelsB = Range_orbB.second;

  XTP_LOG(logDEBUG, *_pLog)
      << "Levels:Basis A[" << levelsA << ":" << basisA << "]"
      << " B[" << levelsB << ":" << basisB << "]" << flush;

  if ((levelsA == 0) || (levelsB == 0)) {
    throw std::runtime_error(
        "No information about number of occupied/unoccupied levels is stored");
  }

  // constructing merged orbitals
  auto MOsA = orbitalsA.MOs().eigenvectors().block(0, Range_orbA.first, basisA,
                                                   Range_orbA.second);
  auto MOsB = orbitalsB.MOs().eigenvectors().block(0, Range_orbB.first, basisB,
                                                   Range_orbB.second);

  XTP_LOG(logDEBUG, *_pLog) << "Calculating overlap matrix for basisset: "
                            << orbitalsAB.getDFTbasisName() << flush;
  Eigen::MatrixXd overlap =
      CalculateOverlapMatrix(orbitalsAB) * orbitalsAB.MOs().eigenvectors();

  XTP_LOG(logDEBUG, *_pLog)
      << "Projecting monomers onto dimer orbitals" << flush;
  Eigen::MatrixXd A_AB = MOsA.transpose() * overlap.topRows(basisA);
  Eigen::MatrixXd B_AB = MOsA.transpose() * overlap.bottomRows(basisB);
  Eigen::VectorXd mag_A = A_AB.rowwise().squaredNorm();
  if (mag_A.any() < 0.95) {
    XTP_LOG(logERROR, *_pLog)
        << "\nWarning: "
        << "Projection of orbitals of monomer A on dimer is insufficient,mag="
        << mag_A.minCoeff() << flush;
  }
  Eigen::VectorXd mag_B = B_AB.rowwise().squaredNorm();
  if (mag_B.any() < 0.95) {
    XTP_LOG(logERROR, *_pLog)
        << "\nWarning: "
        << "Projection of orbitals of monomer B on dimer is insufficient,mag="
        << mag_B.minCoeff() << flush;
  }

  Eigen::MatrixXd psi_AxB_dimer_basis(A_AB.rows() + B_AB.rows(), A_AB.cols());
  psi_AxB_dimer_basis.topRows(A_AB.rows()) = A_AB;
  psi_AxB_dimer_basis.bottomRows(B_AB.rows()) = B_AB;

  XTP_LOG(logDEBUG, *_pLog)
      << "Projecting the Fock matrix onto the dimer basis" << flush;
  Eigen::MatrixXd JAB_dimer = psi_AxB_dimer_basis *
                              orbitalsAB.MOs().eigenvalues().asDiagonal() *
                              psi_AxB_dimer_basis.transpose();
  XTP_LOG(logDEBUG, *_pLog) << "Constructing Overlap matrix" << flush;
  Eigen::MatrixXd S_AxB = psi_AxB_dimer_basis * psi_AxB_dimer_basis.transpose();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_AxB);
  Eigen::MatrixXd Sm1 = es.operatorInverseSqrt();
  XTP_LOG(logDEBUG, *_pLog) << "Smallest eigenvalue of overlap matrix is "
                            << es.eigenvalues()(0) << flush;
  JAB = Sm1 * JAB_dimer * Sm1;
  XTP_LOG(logDEBUG, *_pLog) << "Done with electronic couplings" << flush;
}

}  // namespace xtp
}  // namespace votca
