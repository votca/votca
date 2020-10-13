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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/fourcenter.h"
#include "votca/xtp/aobasis.h"
#include "votca/xtp/make_libint_work.h"
#include <libint2.hpp>

namespace votca {
namespace xtp {

void FCMatrix::Fill_4c_small_molecule(const AOBasis& dftbasis) {

  Index dftBasisSize = dftbasis.AOBasisSize();
  Index vectorSize = (dftBasisSize * (dftBasisSize + 1)) / 2;

  try {
    _4c_vector = Eigen::VectorXd::Zero((vectorSize * (vectorSize + 1)) / 2);
  } catch (std::bad_alloc& ba) {
    throw std::runtime_error(
        "Basisset too large for 4c calculation. Not enough RAM.");
  }
  Index shellsize = dftbasis.getNumofShells();

  // Setup LibInt
  std::vector<libint2::Shell> libintShells = dftbasis.GenerateLibintBasis();

  std::vector<Index> startpos = dftbasis.getMapToBasisFunctions();
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] =
      libint2::Engine(libint2::Operator::coulomb, dftbasis.getMaxNprim(),
                      static_cast<int>(dftbasis.getMaxL()), 0);
  engines[0].set(libint2::BraKet::xx_xx);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

#pragma omp parallel for schedule(dynamic)
  for (Index i = 0; i < shellsize; ++i) {
    libint2::Engine& engine = engines[OPENMP::getThreadId()];
    const libint2::Shell& sh3 = libintShells[i];
    Index start_3 = startpos[i];
    Index NumFunc_3 = sh3.size();

    for (Index j = i; j < shellsize; ++j) {
      const libint2::Shell& sh4 = libintShells[j];
      Index start_4 = startpos[j];
      Index NumFunc_4 = sh4.size();

      for (Index k = i; k < shellsize; ++k) {
        const libint2::Shell& sh1 = libintShells[k];
        Index start_1 = startpos[k];
        Index NumFunc_1 = sh1.size();

        for (Index l = k; l < shellsize; ++l) {
          const libint2::Shell& sh2 = libintShells[l];
          Index start_2 = startpos[l];
          Index NumFunc_2 = sh2.size();

          Eigen::Tensor<double, 4> block;
          bool nonzero =
              FillFourCenterRepBlock(block, engine, sh1, sh2, sh3, sh4);

          if (nonzero) {

            for (Index i_3 = 0; i_3 < NumFunc_3; i_3++) {
              Index ind_3 = start_3 + i_3;
              Index sum_ind_3 = (ind_3 * (ind_3 + 1)) / 2;
              for (Index i_4 = 0; i_4 < NumFunc_4; i_4++) {
                Index ind_4 = start_4 + i_4;
                if (ind_3 > ind_4) {
                  continue;
                }
                Index index_34 = dftBasisSize * ind_3 - sum_ind_3 + ind_4;
                Index index_34_12_a =
                    vectorSize * index_34 - (index_34 * (index_34 + 1)) / 2;
                for (Index i_1 = 0; i_1 < NumFunc_1; i_1++) {
                  Index ind_1 = start_1 + i_1;
                  Index sum_ind_1 = (ind_1 * (ind_1 + 1)) / 2;
                  for (Index i_2 = 0; i_2 < NumFunc_2; i_2++) {
                    Index ind_2 = start_2 + i_2;
                    if (ind_1 > ind_2) {
                      continue;
                    }
                    Index index_12 = dftBasisSize * ind_1 - sum_ind_1 + ind_2;
                    if (index_34 > index_12) {
                      continue;
                    }
                    _4c_vector(index_34_12_a + index_12) =
                        block(i_1, i_2, i_3, i_4);

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

bool FCMatrix::FillFourCenterRepBlock(Eigen::Tensor<double, 4>& block,
                                      libint2::Engine& engine,
                                      const libint2::Shell& shell_1,
                                      const libint2::Shell& shell_2,
                                      const libint2::Shell& shell_3,
                                      const libint2::Shell& shell_4) const {

  const libint2::Engine::target_ptr_vec& buf = engine.results();

  engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
      shell_1, shell_2, shell_3, shell_4);

  Eigen::TensorMap<Eigen::Tensor<const double, 4, Eigen::RowMajor> const>
      result(buf[0], shell_1.size(), shell_2.size(), shell_3.size(),
             shell_4.size());
  std::array<int, 4> shuffle{3, 2, 1, 0};
  // this just turns it into column major order
  block = result.swap_layout().shuffle(shuffle);

  return true;
}

}  // namespace xtp
}  // namespace votca
