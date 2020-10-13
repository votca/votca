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
                                      const AOShell& shell_1,
                                      const AOShell& shell_2,
                                      const AOShell& shell_3,
                                      const AOShell& shell_4) const {

  const libint2::Engine::target_ptr_vec& buf = engine.results();

  libint2::Shell shell1 = shell_1.LibintShell();
  libint2::Shell shell2 = shell_2.LibintShell();
  libint2::Shell shell3 = shell_3.LibintShell();
  libint2::Shell shell4 = shell_4.LibintShell();

  engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
      shell1, shell2, shell3, shell4);

  Eigen::TensorMap<Eigen::Tensor<const double, 4, Eigen::RowMajor> const>
      result(buf[0], shell1.size(), shell2.size(), shell3.size(),
             shell4.size());

  for (size_t index1 = 0; index1 < shell1.size(); index1++) {
    for (size_t index2 = 0; index2 < shell2.size(); index2++) {
      for (size_t index3 = 0; index3 < shell3.size(); index3++) {
        for (size_t index4 = 0; index4 < shell4.size(); index4++) {
          // ao3c is col major result is row-major so this is ideal
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
