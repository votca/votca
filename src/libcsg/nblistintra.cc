/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

// Local VOTCA includes
#include "votca/csg/nblistintra.h"
#include "votca/csg/topology.h"
#include "votca/tools/NDimVector.h"

namespace votca {
namespace csg {

using namespace std;

void NBListIntra::Generate(BeadList &list1, BeadList &list2, bool do_exclusions) {
  BeadList::iterator iter1;
  BeadList::iterator iter2;
  do_exclusions_ = do_exclusions;

  if (list1.empty()) {
    return;
  }
  if (list2.empty()) {
    return;
  }

  assert(&(list1.getTopology()) == &(list2.getTopology()));
  const Topology &top = list1.getTopology();

  for (iter1 = list1.begin(); iter1 != list1.end(); ++iter1) {

    if (&list1 == &list2) {
      iter2 = iter1;
      ++iter2;
    } else {
      iter2 = list2.begin();
    }

    if (iter2 == list2.end()) {
      continue;
    }
    if (*iter1 == *iter2) {
      continue;
    }

    for (; iter2 != list2.end(); ++iter2) {
      // skip if the two beads are not from the same molecule
      if ((*iter1)->getMoleculeId() != (*iter2)->getMoleculeId()) {
        continue;
      }
      Eigen::Vector3d u = (*iter1)->getPos();
      Eigen::Vector3d v = (*iter2)->getPos();

      Eigen::Vector3d r = top.BCShortestConnection(u, v);
      double d = r.norm();
      if (d < cutoff_) {
        if (do_exclusions_) {
          if (top.getExclusions().IsExcluded(*iter1, *iter2)) {
            continue;
          }
        }
        if ((*match_function_)(*iter1, *iter2, r, d)) {
          if (!FindPair(*iter1, *iter2)) {
            AddPair(pair_creator_(*iter1, *iter2, r));
          }
        }
      }
    }
  }
}

}  // namespace csg
}  // namespace votca
