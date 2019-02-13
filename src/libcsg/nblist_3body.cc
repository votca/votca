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

#include <iostream>
#include <votca/csg/nblist_3body.h>

namespace votca {
namespace csg {

NBList_3Body::NBList_3Body() : _do_exclusions(false), _match_function(0) {
  setTripleType<BeadTriple>();
  SetMatchFunction(NBList_3Body::match_always);
}

NBList_3Body::~NBList_3Body() {}

void NBList_3Body::Generate(BeadList &list1, BeadList &list2, BeadList &list3,
                            bool do_exclusions) {
  BeadList::iterator iter1;
  BeadList::iterator iter2;
  BeadList::iterator iter3;
  _do_exclusions = do_exclusions;

  if (list1.empty())
    return;
  if (list2.empty())
    return;
  if (list3.empty())
    return;

  // check if all bead lists "have" the same topology
  assert(list1.getTopology() == list2.getTopology());
  assert(list1.getTopology() == list3.getTopology());
  assert(list2.getTopology() == list3.getTopology());
  Topology *top = list1.getTopology();

  // builds neighbor lists for all cases, where all list are of different bead
  // typess (list1 neq list2 neq list3), list2 and list3 are of the same type
  // (list1 neq (list2 = list3)) or all three lists are the same
  // (list1=list2=list3)!
  for (iter1 = list1.begin(); iter1 != list1.end(); ++iter1) {

    // in all cases iterate over the full bead lists list1 and list2
    // if both lists are the same (the 3 bead types are the same),
    // still necessary as both permutations (1,2,*) (2,1,*) included in neighbor
    // list! (i,jneqi,k>j)
    iter2 = list2.begin();

    for (; iter2 != list2.end(); ++iter2) {

      // do not include the same beads twice in one triple!
      if (*iter1 == *iter2) {
        continue;
      }

      // if list2 and list3 are the same, set iter3 to "iter2+1" (i,jneqi,k>j)
      if (&list2 == &list3) {
        iter3 = iter2;
        ++iter3;
      } else {
        iter3 = list3.begin();
      }

      for (; iter3 != list3.end(); ++iter3) {

        // do not include the same beads twice in one triple!
        if (*iter1 == *iter3) {
          continue;
        }
        if (*iter2 == *iter3) {
          continue;
        }

        vec u = (*iter1)->getPos();
        vec v = (*iter2)->getPos();
        vec z = (*iter3)->getPos();

        vec r12 = top->BCShortestConnection(u, v);
        vec r13 = top->BCShortestConnection(u, z);
        vec r23 = top->BCShortestConnection(v, z);
        double d12 = abs(r12);
        double d13 = abs(r13);
        double d23 = abs(r23);
        // to do: at the moment use only one cutoff value
        // to do: so far only check the distance between bead 1 (central bead)
        // and bead2 and bead 3
        if ((d12 < _cutoff) && (d13 < _cutoff)) {
          /// experimental: at the moment exclude interaction as soon as one of
          /// the three pairs (1,2) (1,3) (2,3) is excluded!
          if (_do_exclusions)
            if ((top->getExclusions().IsExcluded(*iter1, *iter2)) ||
                (top->getExclusions().IsExcluded(*iter1, *iter3)) ||
                (top->getExclusions().IsExcluded(*iter2, *iter3))) {
              continue;
            }
          if ((*_match_function)(*iter1, *iter2, *iter3, r12, r13, r23, d12,
                                 d13, d23))
            if (!FindTriple(*iter1, *iter2, *iter3))
              AddTriple(_triple_creator(*iter1, *iter2, *iter3, r12, r13, r23));
        }
      }
    }
  }
}

} // namespace csg
} // namespace votca
