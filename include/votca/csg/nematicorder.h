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

#ifndef _VOTCA_CSG_NEMATICORDER_H
#define _VOTCA_CSG_NEMATICORDER_H

#include "topology.h"
#include <votca/tools/eigen.h>

namespace votca {
namespace csg {

class NematicOrder {
 public:
  NematicOrder() = default;
  ~NematicOrder() = default;

  void Process(Topology &top, const std::string &filter = "*");

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> &NematicU() {
    return _nemat_u;
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> &NematicV() {
    return _nemat_v;
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> &NematicW() {
    return _nemat_w;
  }

 private:
  Eigen::Matrix3d _mu, _mv, _mw;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> _nemat_u, _nemat_v, _nemat_w;
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_NEMATICORDER_H */
