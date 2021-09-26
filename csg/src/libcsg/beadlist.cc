/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

// VOTCA includes
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/csg/beadlist.h"
#include "votca/csg/topology.h"

namespace votca {
namespace csg {

using namespace std;

Index BeadList::Generate(Topology &top, const string &select) {
  topology_ = &top;
  bool selectByName = false;
  string pSelect;  // parsed selection string

  if (select.substr(0, 5) == "name:") {
    // select according to bead name instead of type
    pSelect = select.substr(5);
    selectByName = true;
  } else {
    pSelect = select;
  }

  for (auto &bead : top.Beads()) {
    if (!selectByName) {
      if (tools::wildcmp(pSelect, bead.getType())) {
        beads_.push_back(&bead);
      }
    } else {
      if (tools::wildcmp(pSelect, bead.getName())) {
        beads_.push_back(&bead);
      }
    }
  }
  return size();
}

Index BeadList::GenerateInSphericalSubvolume(Topology &top,
                                             const string &select,
                                             Eigen::Vector3d ref,
                                             double radius) {
  topology_ = &top;
  bool selectByName = false;
  string pSelect;  // parsed selection string

  if (select.substr(0, 5) == "name:") {
    // select according to bead name instead of type
    pSelect = select.substr(5);
    selectByName = true;
  } else {
    pSelect = select;
  }

  for (auto &bead : top.Beads()) {
    if (topology_->BCShortestConnection(ref, bead.getPos()).norm() > radius) {
      continue;
    }
    if (!selectByName) {
      if (tools::wildcmp(pSelect, bead.getType())) {
        beads_.push_back(&bead);
      }
    } else {
      if (tools::wildcmp(pSelect, bead.getName())) {
        beads_.push_back(&bead);
      }
    }
  }
  return size();
}

}  // namespace csg
}  // namespace votca
