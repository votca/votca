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
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/dftcoupling.h"

namespace votca {
namespace xtp {

using boost::format;
using std::flush;

void DFTcoupling::Initialize(tools::Property& options) {

  degeneracy_ = options.get("degeneracy").as<double>();
  degeneracy_ *= tools::conv::ev2hrt;
  numberofstatesA_ = options.get("levA").as<Index>();
  numberofstatesB_ = options.get("levB").as<Index>();
}

void DFTcoupling::WriteToProperty(tools::Property& type_summary,
                                  const Orbitals& orbitalsA,
                                  const Orbitals& orbitalsB, Index a,
                                  Index b) const {
  double J = getCouplingElement(a, b, orbitalsA, orbitalsB);
  tools::Property& coupling = type_summary.add("coupling", "");
  coupling.setAttribute("levelA", a);
  coupling.setAttribute("levelB", b);
  coupling.setAttribute("j", (format("%1$1.6e") % J).str());
}

void DFTcoupling::Addoutput(tools::Property& type_summary,
                            const Orbitals& orbitalsA,
                            const Orbitals& orbitalsB) const {
  tools::Property& dftcoupling = type_summary.add(Identify(), "");
  dftcoupling.setAttribute("homoA", orbitalsA.getHomo());
  dftcoupling.setAttribute("homoB", orbitalsB.getHomo());
  tools::Property& hole_summary = dftcoupling.add("hole", "");
  // hole hole
  for (Index a = Range_orbA.first; a <= orbitalsA.getHomo(); ++a) {
    for (Index b = Range_orbB.first; b <= orbitalsB.getHomo(); ++b) {
      WriteToProperty(hole_summary, orbitalsA, orbitalsB, a, b);
    }
  }
  tools::Property& electron_summary = dftcoupling.add("electron", "");
  // electron-electron
  for (Index a = orbitalsA.getLumo();
       a <= Range_orbA.first + Range_orbA.second - 1; ++a) {
    for (Index b = orbitalsB.getLumo();
         b <= Range_orbB.first + Range_orbB.second - 1; ++b) {
      WriteToProperty(electron_summary, orbitalsA, orbitalsB, a, b);
    }
  }
  return;
}

std::pair<Index, Index> DFTcoupling::DetermineRangeOfStates(
    const Orbitals& orbital, Index numberofstates) const {
  const Eigen::VectorXd& MOEnergies = orbital.MOs().eigenvalues();
  if (std::abs(MOEnergies(orbital.getHomo()) - MOEnergies(orbital.getLumo())) <
      degeneracy_) {
    throw std::runtime_error(
        "Homo Lumo Gap is smaller than degeneracy. "
        "Either your degeneracy is too large or your Homo and Lumo are "
        "degenerate");
  }

  Index minimal = orbital.getHomo() - numberofstates + 1;
  Index maximal = orbital.getLumo() + numberofstates - 1;

  std::vector<Index> deg_min = orbital.CheckDegeneracy(minimal, degeneracy_);
  minimal = *std::min_element(deg_min.begin(), deg_min.end());

  std::vector<Index> deg_max = orbital.CheckDegeneracy(maximal, degeneracy_);
  maximal = *std::max_element(deg_max.begin(), deg_max.end());

  std::pair<Index, Index> result;
  result.first = minimal;                 // start
  result.second = maximal - minimal + 1;  // size

  return result;
}

double DFTcoupling::getCouplingElement(Index levelA, Index levelB,
                                       const Orbitals& orbitalsA,
                                       const Orbitals& orbitalsB) const {

  Index levelsA = Range_orbA.second;
  if (degeneracy_ != 0) {
    std::vector<Index> list_levelsA =
        orbitalsA.CheckDegeneracy(levelA, degeneracy_);
    std::vector<Index> list_levelsB =
        orbitalsB.CheckDegeneracy(levelB, degeneracy_);

    double JAB_sq = 0;

    for (Index iA : list_levelsA) {
      Index indexA = iA - Range_orbA.first;
      for (Index iB : list_levelsB) {
        Index indexB = iB - Range_orbB.first + levelsA;
        double JAB_one_level = JAB(indexA, indexB);
        JAB_sq += JAB_one_level * JAB_one_level;
      }
    }
    return std::sqrt(JAB_sq /
                     double(list_levelsA.size() * list_levelsB.size())) *
           tools::conv::hrt2ev;
  } else {
    Index indexA = levelA - Range_orbA.first;
    Index indexB = levelB - Range_orbB.first + levelsA;
    return JAB(indexA, indexB) * tools::conv::hrt2ev;
  }
}

/**
 * \brief evaluates electronic couplings
 * @param  orbitalsA molecular orbitals of molecule A
 * @param  orbitalsB molecular orbitals of molecule B
 * @param  orbitalsAB molecular orbitals of the dimer AB
 */
void DFTcoupling::CalculateCouplings(const Orbitals& orbitalsA,
                                     const Orbitals& orbitalsB,
                                     const Orbitals& orbitalsAB) {

  XTP_LOG(Log::error, *pLog_) << "Calculating electronic couplings" << flush;

  CheckAtomCoordinates(orbitalsA, orbitalsB, orbitalsAB);

  // constructing the direct product orbA x orbB
  Index basisA = orbitalsA.getBasisSetSize();
  Index basisB = orbitalsB.getBasisSetSize();

  if ((basisA == 0) || (basisB == 0)) {
    throw std::runtime_error("Basis set size is not stored in monomers");
  }

  Range_orbA = DetermineRangeOfStates(orbitalsA, numberofstatesA_);
  Range_orbB = DetermineRangeOfStates(orbitalsB, numberofstatesB_);

  Index levelsA = Range_orbA.second;
  Index levelsB = Range_orbB.second;

  XTP_LOG(Log::error, *pLog_)
      << "Levels:Basis A[" << levelsA << ":" << basisA << "]" << " B["
      << levelsB << ":" << basisB << "]" << flush;

  if ((levelsA == 0) || (levelsB == 0)) {
    throw std::runtime_error(
        "No information about number of occupied/unoccupied levels is stored");
  }

  // constructing merged orbitals
  auto MOsA = orbitalsA.MOs().eigenvectors().middleCols(Range_orbA.first,
                                                        Range_orbA.second);
  auto MOsB = orbitalsB.MOs().eigenvectors().middleCols(Range_orbB.first,
                                                        Range_orbB.second);

  XTP_LOG(Log::info, *pLog_) << "Calculating overlap matrix for basisset: "
                             << orbitalsAB.getDFTbasisName() << flush;

  Eigen::MatrixXd overlap =
      CalculateOverlapMatrix(orbitalsAB) * orbitalsAB.MOs().eigenvectors();

  XTP_LOG(Log::info, *pLog_)
      << "Projecting monomers onto dimer orbitals" << flush;
  Eigen::MatrixXd A_AB = MOsA.transpose() * overlap.topRows(basisA);
  Eigen::MatrixXd B_AB = MOsB.transpose() * overlap.bottomRows(basisB);
  Eigen::VectorXd mag_A = A_AB.rowwise().squaredNorm();
  if (mag_A.any() < 0.95) {
    XTP_LOG(Log::error, *pLog_)
        << "\nWarning: "
        << "Projection of orbitals of monomer A on dimer is insufficient,mag="
        << mag_A.minCoeff() << flush;
  }
  Eigen::VectorXd mag_B = B_AB.rowwise().squaredNorm();
  if (mag_B.any() < 0.95) {
    XTP_LOG(Log::error, *pLog_)
        << "\nWarning: "
        << "Projection of orbitals of monomer B on dimer is insufficient,mag="
        << mag_B.minCoeff() << flush;
  }

  Eigen::MatrixXd psi_AxB_dimer_basis(A_AB.rows() + B_AB.rows(), A_AB.cols());
  psi_AxB_dimer_basis.topRows(A_AB.rows()) = A_AB;
  psi_AxB_dimer_basis.bottomRows(B_AB.rows()) = B_AB;

  XTP_LOG(Log::info, *pLog_)
      << "Projecting the Fock matrix onto the dimer basis" << flush;
  Eigen::MatrixXd JAB_dimer = psi_AxB_dimer_basis *
                              orbitalsAB.MOs().eigenvalues().asDiagonal() *
                              psi_AxB_dimer_basis.transpose();
  XTP_LOG(Log::info, *pLog_) << "Constructing Overlap matrix" << flush;
  Eigen::MatrixXd S_AxB = psi_AxB_dimer_basis * psi_AxB_dimer_basis.transpose();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_AxB);
  Eigen::MatrixXd Sm1 = es.operatorInverseSqrt();
  XTP_LOG(Log::info, *pLog_) << "Smallest eigenvalue of overlap matrix is "
                             << es.eigenvalues()(0) << flush;
  JAB = Sm1 * JAB_dimer * Sm1;
  XTP_LOG(Log::error, *pLog_) << "Done with electronic couplings" << flush;
}

}  // namespace xtp
}  // namespace votca
