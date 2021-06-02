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

#include "votca/xtp/activedensitymatrix.h"
#include "votca/xtp/aomatrix.h"
#include <votca/tools/eigenio_matrixmarket.h>

namespace votca {
namespace xtp {
void ActiveDensityMatrix::compute_activedensitymatrix(
    Eigen::MatrixXd &new_mo_coeff) {
      QMMolecule mol = orbitals.QMAtoms();
  basis.Load(orbitals.getDFTbasisName());
  aobasis.Fill(basis, mol);
  AOOverlap overlap;
  overlap.Fill(aobasis);
  Eigen::MatrixXd S = overlap.Matrix();
  Index counter = 0;
  std::vector<Index> numfuncpatom = aobasis.getFuncPerAtom();
  Eigen::MatrixXd active_mo_coeff;

  for (Index i = 0; i < new_mo_coeff.cols(); i++) {
    /* calculate <i|P|i> */
    Eigen::MatrixXd multipliedmatrix =
        new_mo_coeff.col(i).transpose() * S * new_mo_coeff.col(i).asDiagonal();
    Eigen::RowVectorXd iP_u_i = multipliedmatrix.colwise().sum();
    Index start = 0;
    for (Index atom_id = 0; atom_id < Index(numfuncpatom.size()); atom_id++) {
      double iPi_x = iP_u_i.segment(start, numfuncpatom[atom_id]).sum();
      if ((std::find(activeatoms.begin(),  activeatoms.end(), atom_id) != activeatoms.end()) && iPi_x > 0.4) {
        active_mo_coeff.conservativeResize(new_mo_coeff.rows(),counter+1);
        active_mo_coeff.col(counter) = new_mo_coeff.col(i);
        counter += 1;
      }
      start += numfuncpatom[atom_id];
    }
  }
  Eigen::MatrixXd Dmat_A = active_mo_coeff * active_mo_coeff.transpose();
}
}  // namespace xtp
}  // namespace votca