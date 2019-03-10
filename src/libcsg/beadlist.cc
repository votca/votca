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

#include <votca/csg/beadlist.h>
#include <votca/csg/topology.h>
#include <votca/tools/tokenizer.h>

namespace votca {
namespace csg {

using namespace std;

int BeadList::Generate(Topology &top, const string &select) {
  BeadContainer::iterator iter;
  _topology           = &top;
  bool   selectByName = false;
  string pSelect;  // parsed selection string

  if (select.substr(0, 5) == "name:") {
    // select according to bead name instead of type
    pSelect      = select.substr(5);
    selectByName = true;
  } else {
    pSelect = select;
  }

  for (iter = top.Beads().begin(); iter != top.Beads().end(); ++iter) {
    if (!selectByName) {
      if (wildcmp(pSelect.c_str(), (*iter)->getType().c_str())) {
        _beads.push_back(*iter);
      }
    } else {
      if (wildcmp(pSelect.c_str(), (*iter)->getName().c_str())) {
        _beads.push_back(*iter);
      }
    }
  }
  return size();
}

int BeadList::GenerateInSphericalSubvolume(Topology &top, const string &select,
                                           Eigen::Vector3d ref, double radius) {
  BeadContainer::iterator iter;
  _topology           = &top;
  bool   selectByName = false;
  string pSelect;  // parsed selection string

  if (select.substr(0, 5) == "name:") {
    // select according to bead name instead of type
    pSelect      = select.substr(5);
    selectByName = true;
  } else {
    pSelect = select;
  }

  for (iter = top.Beads().begin(); iter != top.Beads().end(); ++iter) {
    if (_topology->BCShortestConnection(ref, (*iter)->getPos()).norm() > radius)
      continue;
    if (!selectByName) {
      if (wildcmp(pSelect.c_str(), (*iter)->getType().c_str())) {
        _beads.push_back(*iter);
      }
    } else {
      if (wildcmp(pSelect.c_str(), (*iter)->getName().c_str())) {
        _beads.push_back(*iter);
      }
    }
  }
  return size();
}

}  // namespace csg
}  // namespace votca
