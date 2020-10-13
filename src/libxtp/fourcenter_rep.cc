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
 * distributed under the License is distributed on an "AS iS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/aotransform.h"
#include "votca/xtp/fourcenter.h"
#include <libint2.hpp>

namespace votca {
namespace xtp {

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

  for (size_t index1 = 0; index1 < shell_1.size(); index1++) {
    for (size_t index2 = 0; index2 < shell_2.size(); index2++) {
      for (size_t index3 = 0; index3 < shell_3.size(); index3++) {
        for (size_t index4 = 0; index4 < shell_4.size(); index4++) {
          block(index1, index2, index3, index4) =
              result(index1, index2, index3, index4);
        }
      }
    }
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
