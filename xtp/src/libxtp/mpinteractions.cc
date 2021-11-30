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
 */

#include <votca/xtp/mpinteractions.h>

namespace votca {
namespace xtp {

MPInteractions::MPInteractions(double alpha, double thole_damping) {
  alpha_ = alpha;
  thole = thole_damping;
  thole2 = thole * thole;
  thole3 = thole * thole2;
}

Eigen::Vector3d MPInteractions::r_fieldAtBy(const Multipole& mp,
                                            const Eigen::Vector3d& dr) {
  InteractionTensor<Screening::erfc, 3> T(alpha_);
  T.computeTensors(dr);
  return -mp.charge * T.rank1() + T.rank2() * mp.dipole; 
}

}  // namespace xtp
}  // namespace votca