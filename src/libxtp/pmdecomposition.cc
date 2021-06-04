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

#include "votca/xtp/pmdecomposition.h"
#include "votca/xtp/aomatrix.h"
#include <votca/tools/eigenio_matrixmarket.h>

namespace votca {
namespace xtp {
void PMDecomposition::computePMD(Orbitals &orbitals) {
  Eigen::MatrixXd mo_coeff = orbitals.MOs().eigenvectors().leftCols(
      orbitals.getNumberOfAlphaElectrons());
  QMMolecule mol = orbitals.QMAtoms();
  basis.Load(orbitals.getDFTbasisName());
  aobasis.Fill(basis, mol);
  AOOverlap overlap;
  overlap.Fill(aobasis);
  Eigen::MatrixXd S = overlap.Matrix();
  double diff_D = 10000;
  Index i = 1;
  while (diff_D > 1e-6) {
    XTP_LOG(Log::error, log) << "Iteration: " << i << std::flush;
    Eigen::MatrixXd uppertriangular = orbitalselections(mo_coeff, S);
    Index maxrow, maxcol;
    diff_D = uppertriangular.maxCoeff(&maxrow, &maxcol);
    XTP_LOG(Log::error, log) << maxrow << " " << maxcol << std::flush;
    XTP_LOG(Log::error, log) << diff_D << std::flush;
    Eigen::MatrixXd max_orbs(mo_coeff.rows(), 2);
    max_orbs << mo_coeff.col(maxrow), mo_coeff.col(maxcol);
    Eigen::MatrixXd new_orbs(mo_coeff.rows(), 2);
    new_orbs = rotatedorbitals(max_orbs, maxrow, maxcol);
    update_maximums(mo_coeff, maxrow, maxcol, new_orbs);

    i += 1;
  }
  orbitals.setPMLocalizedOrbitals(mo_coeff);
}

Eigen::MatrixXd PMDecomposition::rotatedorbitals(Eigen::MatrixXd &maxorbs,
                                                 Index s, Index t) {
  Eigen::MatrixXd neworbitals(maxorbs.rows(), 2);
  Eigen::VectorXd vec1 = maxorbs.col(0);
  Eigen::VectorXd vec2 = maxorbs.col(1);
  double gam =
      0.25 * asin(B(s, t) / sqrt((A(s, t) * A(s, t)) + (B(s, t) * B(s, t))));
  double sin_gamma = std::sin(gam);
  double cos_gamma = std::cos(gam);
  Eigen::VectorXd new_vec1 = (cos_gamma * vec1) + (sin_gamma * vec2);
  Eigen::VectorXd new_vec2 = -1 * (sin_gamma * vec1) + (cos_gamma * vec2);
  neworbitals << new_vec1, new_vec2;
  XTP_LOG(Log::error, log) << sin_gamma << std::flush;
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
    for (Index t = 0; t < m.cols(); t++) {
      if (t > s) {
        Eigen::VectorXd req_vec1 = m.col(s);
        Eigen::VectorXd req_vec2 = m.col(t);
        Eigen::MatrixXd a = S * req_vec1.asDiagonal();
        Eigen::MatrixXd b = S * req_vec2.asDiagonal();
        Eigen::MatrixXd c = req_vec1.transpose() * a;   // vec1.S.vec1
        Eigen::MatrixXd d = req_vec1.transpose() * b;   // term1 of eq 31
        Eigen::MatrixXd e = req_vec1.asDiagonal() * b;  // term2 of eq 31
        Eigen::MatrixXd f = req_vec2.transpose() * b;   // vec2.S.vec2
        Eigen::RowVectorXd sps = c.colwise().sum();
        Eigen::RowVectorXd tpt = f.colwise().sum();
        Eigen::RowVectorXd spt =
            0.5 * (d.colwise().sum() + e.rowwise().sum().transpose());
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
  }
  return zeromatrix;
}

void PMDecomposition::update_maximums(Eigen::MatrixXd &m, Index col1,
                                      Index col2, Eigen::MatrixXd &new_orbs) {
  m.col(col1) = new_orbs.col(0);
  m.col(col2) = new_orbs.col(1);
}

}  // namespace xtp
}  // namespace votca