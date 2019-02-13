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
vec TriclinicBox::BCShortestConnection(const vec &r_i, const vec &r_j) const {
  vec r_tp, r_dp, r_sp, r_ij;
  vec a = _box.getCol(0);
  vec b = _box.getCol(1);
  vec c = _box.getCol(2);
  r_tp = r_j - r_i;
  r_dp = r_tp - c * round(r_tp.getZ() / c.getZ());
  r_sp = r_dp - b * round(r_dp.getY() / b.getY());
  r_ij = r_sp - a * round(r_sp.getX() / a.getX());
  return r_ij;
}

} // namespace csg
} // namespace votca
