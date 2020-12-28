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

#include <votca/csg/orthorhombicbox.h>

namespace votca {
namespace csg {

Eigen::Vector3d OrthorhombicBox::BCShortestConnection(
    const Eigen::Vector3d &r_i, const Eigen::Vector3d &r_j) const {
  const Eigen::Array3d box = _box.diagonal();
  const Eigen::Array3d r_ij = r_j - r_i;
  return (r_ij - box * (r_ij / box).round()).matrix();
}

}  // namespace csg
}  // namespace votca
