/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/csg/triclinicbox.h>

namespace votca {
namespace csg {

/*
 This way of determining the shortest distance is only working for a set of
 triclinic boxes, in particular
 a_y = a_z = b_z = 0
 a_x > 0, b_y > 0, c_z > 0
 b_x < 0.5 a_x, c_x < 0.5 a_x, c_y < 0.5 b_y
 GROMACS checks if these conditions are satisfied.
 If a simple search algorithm is used to determine if a particle
 is a within cutoff r_c, make sure that r_c < 0.5 min(a_x, b_y, c_z)
 */
Eigen::Vector3d TriclinicBox::BCShortestConnection(
    const Eigen::Vector3d &r_i, const Eigen::Vector3d &r_j) const {
  Eigen::Vector3d r_tp, r_dp, r_sp, r_ij;
  Eigen::Vector3d a = _box.col(0);
  Eigen::Vector3d b = _box.col(1);
  Eigen::Vector3d c = _box.col(2);
  r_tp              = r_j - r_i;
  r_dp              = r_tp - c * round(r_tp.z() / c.z());
  r_sp              = r_dp - b * round(r_dp.y() / b.y());
  r_ij              = r_sp - a * round(r_sp.x() / a.x());
  return r_ij;
}

}  // namespace csg
}  // namespace votca
