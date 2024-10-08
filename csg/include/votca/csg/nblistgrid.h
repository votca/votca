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

#ifndef VOTCA_CSG_NBLISTGRID_H
#define VOTCA_CSG_NBLISTGRID_H

// Standard includes
#include <vector>

// VOTCA includes
#include <votca/tools/eigen.h>

// Local VOTCA includes
#include "nblist.h"
#include "votca/tools/NDimVector.h"

namespace votca {
namespace csg {

class NBListGrid : public NBList {
 public:
  void Generate(BeadList &list1, BeadList &list2,
                bool do_exclusions = true) override;
  void Generate(BeadList &list, bool do_exclusions = true) override;

 protected:
  struct cell_t {
    BeadList beads_;
    std::vector<cell_t *> neighbours_;
  };

  Eigen::Vector3d box_a_, box_b_, box_c_;
  Eigen::Vector3d norm_a_, norm_b_, norm_c_;
  Index box_Na_, box_Nb_, box_Nc_;

  tools::NDimVector<cell_t, 3> grid_;

  void InitializeGrid(const Eigen::Matrix3d &box);

  cell_t &getCell(const Eigen::Vector3d &r);
  cell_t &getCell(const Index &a, const Index &b, const Index &c);

  void TestBead(const Topology &top, cell_t &cell, Bead *bead);
  void TestCell(const Topology &top, cell_t &cell, Bead *bead);
};

}  // namespace csg
}  // namespace votca

#endif /*  VOTCA_CSG_NBLISTGRID_H */
