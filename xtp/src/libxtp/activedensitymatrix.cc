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

std::array<Eigen::MatrixXd, 2> ActiveDensityMatrix::compute_Dmat_A() {
  Eigen::MatrixXd localized_mo_coeff = orbitals_.getPMLocalizedOrbital();
  std::array<Eigen::MatrixXd, 2> two_distinct_regions =
      activedensitymatrix(localized_mo_coeff);
  return two_distinct_regions;
}

bool ActiveDensityMatrix::checker_function(Eigen::MatrixXd mat,
                                           Eigen::VectorXd vec) {
  bool result = true;
  for (Index i = 0; i < mat.cols(); i++) {
    if (mat.col(i) == vec) {
      result = false;
      break;
    }
  }
  return result;
}

std::array<Eigen::MatrixXd, 2> ActiveDensityMatrix::activedensitymatrix(
    Eigen::MatrixXd &localized_mo_coeff) {
  AOBasis aobasis;
  aobasis = orbitals_.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(aobasis);
  Index counter = 0;
  std::vector<Index> numfuncpatom = aobasis.getFuncPerAtom();
  Eigen::MatrixXd active_mo_coeff;

  for (Index LocMoCoeff_col_i = 0; LocMoCoeff_col_i < localized_mo_coeff.cols();
       LocMoCoeff_col_i++) {
    /* calculate <i|P|i> */
    Eigen::MatrixXd orbital_wise_population =
        localized_mo_coeff.col(LocMoCoeff_col_i).transpose() *
        overlap.Matrix() *
        localized_mo_coeff.col(LocMoCoeff_col_i).asDiagonal();
    Eigen::RowVectorXd iP_u_i = orbital_wise_population.colwise().sum();
    Index start = 0;
    for (Index atom_id = 0; atom_id < Index(numfuncpatom.size()); atom_id++) {
      double iPi_x = iP_u_i.segment(start, numfuncpatom[atom_id]).sum();
      if ((std::find(activeatoms_.begin(), activeatoms_.end(), atom_id) !=
           activeatoms_.end()) &&
          iPi_x > 0.40) {
        std::cout << "Contri of orbital " << LocMoCoeff_col_i << "on atom "
                  << atom_id << "is : " << iPi_x << std::endl;
        if (checker_function(active_mo_coeff,
                             localized_mo_coeff.col(LocMoCoeff_col_i))) {
          active_mo_coeff.conservativeResize(localized_mo_coeff.rows(),
                                             counter + 1);
          active_mo_coeff.col(counter) =
              localized_mo_coeff.col(LocMoCoeff_col_i);
          counter++;
        }
      }
      start += numfuncpatom[atom_id];
    }
  }
  std::array<Eigen::MatrixXd, 2> result;
  std::cout << active_mo_coeff.cols() << std::endl;
  result[0] = 2 * active_mo_coeff * active_mo_coeff.transpose();
  result[1] = localized_mo_coeff;
  std::cout << "No. of electrons: "
            << result[0].cwiseProduct(overlap.Matrix()).sum();
  return result;
}
}  // namespace xtp
}  // namespace votca