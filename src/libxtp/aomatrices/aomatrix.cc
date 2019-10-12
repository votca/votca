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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <vector>
#include <votca/xtp/aomatrix.h>

namespace votca {
namespace xtp {

void AOMatrix::Fill(const AOBasis& aobasis) {
  _aomatrix =
      Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
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
      Eigen::Block<Eigen::MatrixXd> block = _aomatrix.block(
          row_start, col_start, shell_row.getNumFunc(), shell_col.getNumFunc());
      // Fill block
      FillBlock(block, shell_row, shell_col);
    }
  }
  // Fill whole matrix by copying
  _aomatrix = _aomatrix.template selfadjointView<Eigen::Lower>();
  return;
}

}  // namespace xtp
}  // namespace votca
