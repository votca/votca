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
 * Reference- A Simple, Exact Density-Functional-Theory Embedding Scheme
 *Frederick R. Manby, Martina Stella, Jason D. Goodpaster, and Thomas F. Miller
 *Journal of Chemical Theory and Computation 2012 8 (8), 2564-2568
 *DOI: 10.1021/ct300544e
 */

#include "votca/xtp/activedensitymatrix.h"
#include "votca/xtp/aomatrix.h"
namespace votca {
namespace xtp {

Eigen::MatrixXd ActiveDensityMatrix::compute_Dmat_A() {
  Eigen::MatrixXd localized_mo_coeff = orbitals_.getPMLocalizedOrbitals();
  return activedensitymatrix(localized_mo_coeff);
}
Eigen::MatrixXd ActiveDensityMatrix::activedensitymatrix(
    Eigen::MatrixXd &localized_mo_coeff) {
  AOBasis aobasis;
  aobasis = orbitals_.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(aobasis);
  Index counter = 0;
  std::vector<Index> numfuncpatom = aobasis.getFuncPerAtom();
  Eigen::MatrixXd active_mo_coeff;

  for (Index i = 0; i < localized_mo_coeff.cols(); i++) {
    /* calculate <i|P|i> */
    Eigen::MatrixXd multipliedmatrix =
        localized_mo_coeff.col(i).transpose() * overlap.Matrix() * localized_mo_coeff.col(i).asDiagonal();
    Eigen::RowVectorXd iP_u_i = multipliedmatrix.colwise().sum();
    Index start = 0;
    for (Index atom_id = 0; atom_id < Index(numfuncpatom.size()); atom_id++) {
      double iPi_x = iP_u_i.segment(start, numfuncpatom[atom_id]).sum();
      if ((std::find(activeatoms_.begin(), activeatoms_.end(), atom_id) !=
           activeatoms_.end()) &&
          iPi_x > 0.4) {
        active_mo_coeff.conservativeResize(localized_mo_coeff.rows(), counter + 1);
        active_mo_coeff.col(counter) = localized_mo_coeff.col(i);
        counter ++;
      }
      start += numfuncpatom[atom_id];
    }
  }
  return active_mo_coeff * active_mo_coeff.transpose();
}
}  // namespace xtp
}  // namespace votca