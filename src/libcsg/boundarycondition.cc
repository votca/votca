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
#include <cassert>
#include <vector>
#include <votca/csg/boundarycondition.h>

namespace votca {
  namespace csg {

    double BoundaryCondition::BoxVolume() const noexcept {
      return std::abs(_box.determinant());
    }

    double BoundaryCondition::getShortestBoxDimension() const {
      assert(getBoxType() != eBoxtype::typeOpen &&
          "Cannot get the shortest dimension of the box because it is open");

      Eigen::Vector3d box_a = box_.col(0);
      Eigen::Vector3d box_b = box_.col(1);
      Eigen::Vector3d box_c = box_.col(2);

      // create plane normals
      Eigen::Vector3d norm_a = box_b.cross(box_c);
      Eigen::Vector3d norm_b = box_c.cross(box_a);
      Eigen::Vector3d norm_c = box_a.cross(box_b);

      norm_a.normalize();
      norm_b.normalize();
      norm_c.normalize();

      double la = box_a.dot(norm_a);
      double lb = box_b.dot(norm_b);
      double lc = box_c.dot(norm_c);

      return std::min(la, std::min(lb, lc));
    }

}  // namespace csg
}  // namespace votca
