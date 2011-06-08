/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

namespace votca { namespace csg {

vec TriclinicBox::BCShortestConnection(const vec &r_i, const vec &r_j) const
{
    vec r_tp, r_dp, r_sp, r_ij;
    vec a = _box.getCol(0); vec b = _box.getCol(1); vec c = _box.getCol(2);
    r_tp = r_j - r_i;
    r_dp = r_tp - c*round(r_tp.getZ()/c.getZ());
    r_sp = r_dp - b*round(r_dp.getY()/b.getY());
    r_ij = r_sp - a*round(r_sp.getX()/a.getX());
    return r_ij;
}

}}
