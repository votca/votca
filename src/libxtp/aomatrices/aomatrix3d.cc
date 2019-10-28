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

#include <votca/xtp/aomatrix3d.h>

namespace votca {
namespace xtp {

void AOMatrix3D::Fill(const AOBasis& aobasis) {
  for (Index i = 0; i < 3; i++) {
    _aomatrix[i] =
        Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
  }
  // loop row
#pragma omp parallel for
  for (Index row = 0; row < aobasis.getNumofShells(); row++) {
    const AOShell& shell_row = aobasis.getShell(row);
    Index row_start = shell_row.getStartIndex();
    // loop column
    for (const AOShell& shell_col : aobasis) {
      // figure out the submatrix
      Index col_start = shell_col.getStartIndex();
      std::vector<Eigen::Block<Eigen::MatrixXd> > submatrix;
      for (Index i = 0; i < 3; i++) {
        Eigen::Block<Eigen::MatrixXd> block =
            _aomatrix[i].block(row_start, col_start, shell_row.getNumFunc(),
                               shell_col.getNumFunc());
        submatrix.push_back(block);
      }
      // Fill block
      FillBlock(submatrix, shell_row, shell_col);
    }
  }
  return;
}

}  // namespace xtp
}  // namespace votca
