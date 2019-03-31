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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/fourcenter.h>

namespace votca {
namespace xtp {

void FCMatrix::Fill_4c_small_molecule(const AOBasis& dftbasis) {
  tensor4d::extent_gen extents;
  int dftBasisSize = dftbasis.AOBasisSize();
  int vectorSize = (dftBasisSize * (dftBasisSize + 1)) / 2;

  try {
    _4c_vector = Eigen::VectorXd::Zero((vectorSize * (vectorSize + 1)) / 2);
  } catch (std::bad_alloc& ba) {
    throw std::runtime_error(
        "Basisset too large for 4c calculation. Not enough RAM.");
  }
  int shellsize = dftbasis.getNumofShells();
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < shellsize; ++i) {

    const AOShell& shell_3 = dftbasis.getShell(i);
    int start_3 = shell_3.getStartIndex();
    int NumFunc_3 = shell_3.getNumFunc();

    for (int j = i; j < shellsize; ++j) {
      const AOShell& shell_4 = dftbasis.getShell(j);
      int start_4 = shell_4.getStartIndex();
      int NumFunc_4 = shell_4.getNumFunc();

      for (int k = i; k < shellsize; ++k) {
        const AOShell& shell_1 = dftbasis.getShell(k);
        int start_1 = shell_1.getStartIndex();
        int NumFunc_1 = shell_1.getNumFunc();

        for (int l = k; l < shellsize; ++l) {
          const AOShell& shell_2 = dftbasis.getShell(l);
          int start_2 = shell_2.getStartIndex();
          int NumFunc_2 = shell_2.getNumFunc();

          tensor4d block(extents[range(0, NumFunc_1)][range(0, NumFunc_2)]
                                [range(0, NumFunc_3)][range(0, NumFunc_4)]);
          for (int i = 0; i < shell_1.getNumFunc(); ++i) {
            for (int j = 0; j < shell_2.getNumFunc(); ++j) {
              for (int k = 0; k < shell_3.getNumFunc(); ++k) {
                for (int l = 0; l < shell_4.getNumFunc(); ++l) {
                  block[i][j][k][l] = 0.0;
                }
              }
            }
          }
          bool nonzero =
              FillFourCenterRepBlock(block, shell_1, shell_2, shell_3, shell_4);

          if (nonzero) {

            for (int i_3 = 0; i_3 < NumFunc_3; i_3++) {
              int ind_3 = start_3 + i_3;
              int sum_ind_3 = (ind_3 * (ind_3 + 1)) / 2;
              for (int i_4 = 0; i_4 < NumFunc_4; i_4++) {
                int ind_4 = start_4 + i_4;
                if (ind_3 > ind_4) continue;
                int index_34 = dftBasisSize * ind_3 - sum_ind_3 + ind_4;
                int index_34_12_a =
                    vectorSize * index_34 - (index_34 * (index_34 + 1)) / 2;
                for (int i_1 = 0; i_1 < NumFunc_1; i_1++) {
                  int ind_1 = start_1 + i_1;
                  int sum_ind_1 = (ind_1 * (ind_1 + 1)) / 2;
                  for (int i_2 = 0; i_2 < NumFunc_2; i_2++) {
                    int ind_2 = start_2 + i_2;
                    if (ind_1 > ind_2) continue;
                    int index_12 = dftBasisSize * ind_1 - sum_ind_1 + ind_2;
                    if (index_34 > index_12) continue;
                    _4c_vector(index_34_12_a + index_12) =
                        block[i_1][i_2][i_3][i_4];

                  }  // i_2
                }    // i_1
              }      // i_4
            }        // i_3

          }  // end if
        }    // DFT shell_2
      }      // DFT shell_1
    }        // DFT shell_4
  }          // DFT shell_3

  return;
}  // FCMatrix_dft::Fill_4c_small_molecule

}  // namespace xtp
}  // namespace votca
