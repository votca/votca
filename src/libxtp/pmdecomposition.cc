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
 * Reference- A fast intrinsic localization procedure applicable for ab initio and semiempirical linear combination of atomic orbital wave functions
 * J. Chem. Phys. 90, 4916 (1989); https://doi.org/10.1063/1.456588
 * János Pipek and Paul G. Mezey
 */

#include "votca/xtp/pmdecomposition.h"
#include "votca/xtp/aomatrix.h"
//#include <votca/tools/eigenio_matrixmarket.h>
#include <limits>

namespace votca {
namespace xtp {
void PMDecomposition::computePMD(Orbitals &orbitals_) {
  Eigen::MatrixXd occ_orbitals = orbitals_.MOs().eigenvectors().leftCols(
      orbitals_.getNumberOfAlphaElectrons());
  QMMolecule mol = orbitals_.QMAtoms();
  basis.Load(orbitals_.getDFTbasisName());
  aobasis.Fill(basis, mol);
  AOOverlap overlap;
  overlap.Fill(aobasis);
  double diff_D = std::numeric_limits<double>::max();
  Index i = 1;
  while (diff_D > 1e-6 && i < 10000) {
    XTP_LOG(Log::error, log_) << "Iteration: " << i << std::flush;
    Eigen::MatrixXd uppertriangular = orbitalselections(occ_orbitals, overlap.Matrix());
    Index maxrow, maxcol;
    diff_D = uppertriangular.maxCoeff(&maxrow, &maxcol);
    XTP_LOG(Log::error, log_) << "Orbitals to be changed: " << maxrow << " " << maxcol << std::flush;
    XTP_LOG(Log::error, log_) << "change in the penalty function: " << diff_D << std::flush;
    Eigen::MatrixX2d max_orbs(occ_orbitals.rows(), 2);
    max_orbs << occ_orbitals.col(maxrow), occ_orbitals.col(maxcol);
    Eigen::MatrixX2d new_orbs(occ_orbitals.rows(), 2);
    new_orbs = rotatedorbitals(max_orbs, maxrow, maxcol);
    update_maximums(occ_orbitals, maxrow, maxcol, new_orbs);

    i += 1;
  }
  orbitals_.setPMLocalizedOrbitals(occ_orbitals);
}

Eigen::MatrixX2d PMDecomposition::rotatedorbitals(Eigen::MatrixX2d &maxorbs,
                                                 Index s, Index t) {  
  Eigen::MatrixX2d neworbitals(maxorbs.rows(), 2);
  Eigen::VectorXd vec1 = maxorbs.col(0);
  Eigen::VectorXd vec2 = maxorbs.col(1);
  double gam = 0.25 * asin(B(s, t) / sqrt((A(s, t) * A(s, t)) + (B(s, t) * B(s, t))));
  Eigen::VectorXd new_vec1 = (std::cos(gam) * vec1) + (std::sin(gam) * vec2);
  Eigen::VectorXd new_vec2 = -1 * (std::sin(gam) * vec1) + (std::cos(gam) * vec2);
  neworbitals << new_vec1, new_vec2;
  XTP_LOG(Log::error, log_) << "Sine of the rotation angle = " << std::sin(gam) << std::flush;
  return neworbitals;
}

// Function to select n(n-1)/2 orbitals and process Ast and Bst
Eigen::MatrixXd PMDecomposition::orbitalselections(Eigen::MatrixXd &m,
                                                   const Eigen::MatrixXd &S) {
  Eigen::MatrixXd req_mat(m.rows(), 2);
  Eigen::MatrixXd zeromatrix = Eigen::MatrixXd::Zero(m.cols(), m.cols());
  A = Eigen::MatrixXd::Zero(m.cols(), m.cols());
  B = Eigen::MatrixXd::Zero(m.cols(), m.cols());
  for (Index s = 0; s < m.cols(); s++) {
    for (Index t = s+1; t < m.cols(); t++) {
        Eigen::RowVectorXd sps = (m.col(s).asDiagonal() * S * m.col(s).asDiagonal()).colwise().sum();   // vec1.S.vec1
        Eigen::MatrixXd spt_split = m.col(s).asDiagonal() * S * m.col(t).asDiagonal();   // terms of eq 31 referenced above
        Eigen::RowVectorXd tpt = (m.col(t).asDiagonal() * S * m.col(t).asDiagonal()).colwise().sum();   // vec2.S.vec2
        Eigen::RowVectorXd spt = 0.5 * (spt_split.colwise().sum() + spt_split.rowwise().sum().transpose());
        std::vector<Index> numfuncpatom = aobasis.getFuncPerAtom();
        Index start = 0;
        double Ast = 0;
        double Bst = 0;
        for (Index atom_id = 0; atom_id < Index(numfuncpatom.size());
             atom_id++) {
          double sps_x = sps.segment(start, numfuncpatom[atom_id]).sum();
          double spt_x = spt.segment(start, numfuncpatom[atom_id]).sum();
          double tpt_x = tpt.segment(start, numfuncpatom[atom_id]).sum();
          Ast += spt_x * spt_x - 0.25 * ((sps_x - tpt_x) * (sps_x - tpt_x));
          Bst += spt_x * (sps_x - tpt_x);
          start += numfuncpatom[atom_id];
        }
        A(s, t) = Ast;
        B(s, t) = Bst;
        double parameter = Ast + sqrt((Ast * Ast) + (Bst * Bst));
        zeromatrix(s, t) = parameter;
    }
  }
  return zeromatrix;
}

void PMDecomposition::update_maximums(Eigen::MatrixXd &m, Index col1,
                                      Index col2, Eigen::MatrixX2d &new_orbs) {
  m.col(col1) = new_orbs.col(0);
  m.col(col2) = new_orbs.col(1);
}
}  // namespace xtp
}  // namespace votca