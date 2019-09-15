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

#include <votca/xtp/aopotential.h>

namespace votca {
namespace xtp {

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> AOPotential<T>::Fill(
    const AOBasis& aobasis) const {
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXcdd;

  MatrixXcdd result =
      MatrixXcdd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
  // AOMatrix is symmetric, restrict explicit calculation of lower triangular
  // matrix
#pragma omp parallel for schedule(guided)
  for (int col = 0; col < aobasis.getNumofShells(); col++) {
    const AOShell& shell_col = aobasis.getShell(col);
    int col_start = shell_col.getStartIndex();
    for (int row = col; row < aobasis.getNumofShells(); row++) {
      const AOShell& shell_row = aobasis.getShell(row);
      int row_start = shell_row.getStartIndex();
      // figure out the submatrix
      Eigen::Block<MatrixXcdd> block = result.block(
          row_start, col_start, shell_row.getNumFunc(), shell_col.getNumFunc());
      // Fill block
      FillBlock(block, shell_row, shell_col);
    }
  }
  // Fill whole matrix by copying
  result.template triangularView<Eigen::StrictlyUpper>() =
      result.template triangularView<Eigen::StrictlyLower>().adjoint();
  return result;
}

template class AOPotential<double>;
template class AOPotential<std::complex<double> >;

}  // namespace xtp
}  // namespace votca
