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
#include "votca/xtp/aomatrix.h"

namespace votca {
namespace xtp {

void AODipole::Fill(const AOBasis& aobasis) {
  auto results = computeOneBodyIntegrals<libint2::Operator::emultipole1,
                                         std::array<libint2::Shell::real_t, 3>>(
      aobasis, _r);

  for (Index i = 0; i < 3; i++) {
    _aomatrix[i] = results[1 + i];  // emultipole1 returns: overlap, x-dipole,
                                    // y-dipole, z-dipole
  }
}
}  // namespace xtp
}  // namespace votca
