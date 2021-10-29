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
 * Reference- A fast intrinsic localization procedure applicable for ab initio
 * and semiempirical linear combination of atomic orbital wave functions J.
 * Chem. Phys. 90, 4916 (1989); https://doi.org/10.1063/1.456588 János Pipek and
 * Paul G. Mezey
 */

#include "votca/xtp/pmdecomposition.h"
#include "votca/xtp/aomatrix.h"
#include <limits>

namespace votca {
namespace xtp {
void PMDecomposition::computePMD(Orbitals &orbitals) {
  Eigen::MatrixXd occupied_orbitals = orbitals.MOs().eigenvectors().leftCols(
      orbitals.getNumberOfAlphaElectrons());
  aobasis = orbitals.getDftBasis();
  AOOverlap overlap;
  overlap.Fill(aobasis);
  double convergence_limit = std::numeric_limits<double>::max();
  Index iteration = 1;
  while (convergence_limit > 1e-6 && iteration < 1000) {
    XTP_LOG(Log::error, log_) << "Iteration: " << iteration << std::flush;
    Eigen::MatrixXd orbital_pair_function_value =
        orbitalselections(occupied_orbitals, overlap.Matrix());
    Index maxrow, maxcol;
    convergence_limit = orbital_pair_function_value.maxCoeff(&maxrow, &maxcol);
    XTP_LOG(Log::error, log_)
        << "Orbitals to be changed: " << maxrow << " " << maxcol << std::flush;
    XTP_LOG(Log::error, log_)
        << "change in the penalty function: " << convergence_limit
        << std::flush;
    Eigen::MatrixX2d max_orbs(occupied_orbitals.rows(), 2);
    max_orbs << occupied_orbitals.col(maxrow), occupied_orbitals.col(maxcol);
    Eigen::MatrixX2d rotated_orbs = rotateorbitals(max_orbs, maxrow, maxcol);
    occupied_orbitals.col(maxrow) = rotated_orbs.col(0);
    occupied_orbitals.col(maxcol) = rotated_orbs.col(1);
    iteration++;
  }
  orbitals.setPMLocalizedOrbital(occupied_orbitals);
}

// Function to rotate the 2 maximum orbitals (s and t)
Eigen::MatrixX2d PMDecomposition::rotateorbitals(
    const Eigen::MatrixX2d &maxorbs, const Index s, const Index t) {
  const double gamma =
      0.25 * asin(B(s, t) / sqrt((A(s, t) * A(s, t)) + (B(s, t) * B(s, t))));
  Eigen::MatrixX2d rotatedorbitals(maxorbs.rows(), 2);
  rotatedorbitals.col(0) =
      (std::cos(gamma) * maxorbs.col(0)) + (std::sin(gamma) * maxorbs.col(1));
  rotatedorbitals.col(1) = -1 * (std::sin(gamma) * maxorbs.col(0)) +
                       (std::cos(gamma) * maxorbs.col(1));
  XTP_LOG(Log::error, log_)
      << "Sine of the rotation angle = " << std::sin(gamma) << std::flush;
  return rotatedorbitals;
}

// Function to select n(n-1)/2 orbitals and process Ast and Bst as described in paper
Eigen::MatrixXd PMDecomposition::orbitalselections(
    Eigen::MatrixXd &occupied_orbitals, const Eigen::MatrixXd &overlap) {
  Eigen::MatrixXd MullikenPop_all_orbitals =
      Eigen::MatrixXd::Zero(occupied_orbitals.cols(), occupied_orbitals.cols());
  A = Eigen::MatrixXd::Zero(occupied_orbitals.cols(), occupied_orbitals.cols());
  B = Eigen::MatrixXd::Zero(occupied_orbitals.cols(), occupied_orbitals.cols());
  for (Index s = 0; s < occupied_orbitals.cols(); s++) {
    for (Index t = s + 1; t < occupied_orbitals.cols(); t++) {
      Eigen::RowVectorXd MullikenPop_orb_S_per_basis =
          (occupied_orbitals.col(s).asDiagonal() * overlap *
           occupied_orbitals.col(s).asDiagonal())
              .colwise()
              .sum(); 
      Eigen::MatrixXd splitwiseMullikenPop_orb_SandT =
          occupied_orbitals.col(s).asDiagonal() * overlap *
          occupied_orbitals.col(t).asDiagonal(); 
      Eigen::RowVectorXd MullikenPop_orb_T_per_basis =
          (occupied_orbitals.col(t).asDiagonal() * overlap *
           occupied_orbitals.col(t).asDiagonal())
              .colwise()
              .sum();
      Eigen::RowVectorXd MullikenPop_orb_SandT_per_basis = 0.5 * (splitwiseMullikenPop_orb_SandT.colwise().sum() +
                                      splitwiseMullikenPop_orb_SandT.rowwise().sum().transpose());
      std::vector<Index> numfuncpatom = aobasis.getFuncPerAtom();
      Index start = 0;
      double Ast = 0;
      double Bst = 0;
      for (Index atom_id = 0; atom_id < Index(numfuncpatom.size()); atom_id++) {
        double MullikenPop_orb_S_per_atom =
            MullikenPop_orb_S_per_basis.segment(start, numfuncpatom[atom_id])
                .sum();
        double MullikenPop_orb_SandT_per_atom = MullikenPop_orb_SandT_per_basis.segment(start, numfuncpatom[atom_id]).sum();
        double MullikenPop_orb_T_per_atom =
            MullikenPop_orb_T_per_basis.segment(start, numfuncpatom[atom_id])
                .sum();
        Ast +=
            MullikenPop_orb_SandT_per_atom * MullikenPop_orb_SandT_per_atom -
            0.25 * ((MullikenPop_orb_S_per_atom - MullikenPop_orb_T_per_atom) *
                    (MullikenPop_orb_S_per_atom - MullikenPop_orb_T_per_atom));
        Bst +=
            MullikenPop_orb_SandT_per_atom * (MullikenPop_orb_S_per_atom - MullikenPop_orb_T_per_atom);
        start += numfuncpatom[atom_id];
      }
      A(s, t) = Ast;
      B(s, t) = Bst;
      MullikenPop_all_orbitals(s, t) = Ast + sqrt((Ast * Ast) + (Bst * Bst));
    }
  }
  return MullikenPop_all_orbitals;
}

}  // namespace xtp
}  // namespace votca