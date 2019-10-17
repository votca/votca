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

#include <votca/csg/nblistgrid_3body.h>
#include <votca/csg/topology.h>

namespace votca {
namespace csg {

using namespace std;

void NBListGrid_3Body::Generate(BeadList &list1, BeadList &list2,
                                BeadList &list3, bool do_exclusions) {
  BeadList::iterator iter;
  _do_exclusions = do_exclusions;
  if (list1.empty()) {
    return;
  }
  if (list2.empty()) {
    return;
  }
  if (list3.empty()) {
    return;
  }

  // check if all bead lists "have" the same topology
  assert(list1.getTopology() == list2.getTopology());
  assert(list1.getTopology() == list3.getTopology());
  assert(list2.getTopology() == list3.getTopology());
  Topology *top = _top = list1.getTopology();

  InitializeGrid(top->getBox());

  // Add all beads of list1 to _beads1
  for (iter = list1.begin(); iter != list1.end(); ++iter) {
    getCell((*iter)->getPos())._beads1.push_back(*iter);
  }

  // Add all beads of list2 to _beads2
  for (iter = list2.begin(); iter != list2.end(); ++iter) {
    getCell((*iter)->getPos())._beads2.push_back(*iter);
  }

  // Add all beads of list2 to _beads3
  for (iter = list3.begin(); iter != list3.end(); ++iter) {
    getCell((*iter)->getPos())._beads3.push_back(*iter);
  }

  // loop over beads of list 1 again to get the correlations
  for (iter = list1.begin(); iter != list1.end(); ++iter) {
    cell_t &cell = getCell((*iter)->getPos());
    TestBead(cell, *iter);
  }
}

void NBListGrid_3Body::Generate(BeadList &list1, BeadList &list2,
                                bool do_exclusions) {
  BeadList::iterator iter;
  _do_exclusions = do_exclusions;
  if (list1.empty()) {
    return;
  }
  if (list2.empty()) {
    return;
  }

  // check if both bead lists "have" the same topology
  assert(list1.getTopology() == list2.getTopology());
  Topology *top = _top = list1.getTopology();

  InitializeGrid(top->getBox());

  // Add all beads of list1 to _beads1
  for (iter = list1.begin(); iter != list1.end(); ++iter) {
    getCell((*iter)->getPos())._beads1.push_back(*iter);
  }

  // Add all beads of list2 to _beads2
  for (iter = list2.begin(); iter != list2.end(); ++iter) {
    getCell((*iter)->getPos())._beads2.push_back(*iter);
  }

  // In this case type2 and type3 are the same
  for (auto &iter : _grid) {
    iter._beads3 = iter._beads2;
  }

  // loop over beads of list 1 again to get the correlations
  for (iter = list1.begin(); iter != list1.end(); ++iter) {
    cell_t &cell = getCell((*iter)->getPos());
    TestBead(cell, *iter);
  }
}

void NBListGrid_3Body::Generate(BeadList &list, bool do_exclusions) {
  BeadList::iterator iter;
  _do_exclusions = do_exclusions;
  if (list.empty()) {
    return;
  }

  Topology *top = _top = list.getTopology();

  InitializeGrid(top->getBox());

  // Add all beads of list to all! bead lists of the cell
  for (iter = list.begin(); iter != list.end(); ++iter) {
    getCell((*iter)->getPos())._beads1.push_back(*iter);
  }

  for (auto &iter : _grid) {
    iter._beads2 = iter._beads1;
    iter._beads3 = iter._beads1;
  }

  // loop over beads again to get the correlations (as all of the same type
  // here)
  for (iter = list.begin(); iter != list.end(); ++iter) {
    cell_t &cell = getCell((*iter)->getPos());
    TestBead(cell, *iter);
  }
}

void NBListGrid_3Body::InitializeGrid(const Eigen::Matrix3d &box) {
  _box_a = box.col(0);
  _box_b = box.col(1);
  _box_c = box.col(2);

  // create plane normals
  _norm_a = _box_b.cross(_box_c);
  _norm_b = _box_c.cross(_box_a);
  _norm_c = _box_a.cross(_box_b);

  _norm_a.normalize();
  _norm_b.normalize();
  _norm_c.normalize();

  double la = _box_a.dot(_norm_a);
  double lb = _box_b.dot(_norm_b);
  double lc = _box_c.dot(_norm_c);

  // calculate grid size, each grid has to be at least size of cut-off
  _box_Na = max((int)(fabs(la / _cutoff)), 1);
  _box_Nb = max((int)(fabs(lb / _cutoff)), 1);
  _box_Nc = max((int)(fabs(lc / _cutoff)), 1);

  _norm_a = _norm_a / la * (double)_box_Na;
  _norm_b = _norm_b / lb * (double)_box_Nb;
  _norm_c = _norm_c / lc * (double)_box_Nc;

  _grid.resize(_box_Na * _box_Nb * _box_Nc);

  int a1, a2, b1, b2, c1, c2;

  a1 = b1 = c1 = -1;
  a2 = b2 = c2 = 1;

  if (_box_Na < 3) {
    a2 = 0;
  }
  if (_box_Nb < 3) {
    b2 = 0;
  }
  if (_box_Nc < 3) {
    c2 = 0;
  }

  if (_box_Na < 2) {
    a1 = 0;
  }
  if (_box_Nb < 2) {
    b1 = 0;
  }
  if (_box_Nc < 2) {
    c1 = 0;
  }

  // wow, setting up the neighbours is an ugly for construct!
  // loop from N..2*N to avoid if and only use %
  for (int a = _box_Na; a < 2 * _box_Na; ++a) {
    for (int b = _box_Nb; b < 2 * _box_Nb; ++b) {
      for (int c = _box_Nc; c < 2 * _box_Nc; ++c) {
        cell_t &cell = getCell(a % _box_Na, b % _box_Nb, c % _box_Nc);
        for (int aa = a + a1; aa <= a + a2; ++aa) {
          for (int bb = b + b1; bb <= b + b2; ++bb) {
            for (int cc = c + c1; cc <= c + c2; ++cc) {
              // test: for 3body algorithm: each cell is a neighbor of its own
              // !!!
              cell._neighbours.push_back(
                  &getCell(aa % _box_Na, bb % _box_Nb, cc % _box_Nc));
            }
          }
        }
      }
    }
  }
}

NBListGrid_3Body::cell_t &NBListGrid_3Body::getCell(const Eigen::Vector3d &r) {
  int a = (int)floor(r.dot(_norm_a));
  int b = (int)floor(r.dot(_norm_b));
  int c = (int)floor(r.dot(_norm_c));

  if (a < 0) {
    a = _box_Na + a % _box_Na;
  }
  a %= _box_Na;

  if (b < 0) {
    b = _box_Nb + b % _box_Nb;
  }
  b %= _box_Nb;

  if (c < 0) {
    c = _box_Nc + c % _box_Nc;
  }
  c %= _box_Nc;

  return getCell(a, b, c);
}

void NBListGrid_3Body::TestBead(NBListGrid_3Body::cell_t &cell, Bead *bead) {
  BeadList::iterator iter2;
  BeadList::iterator iter3;
  Eigen::Vector3d u = bead->getPos();

  // loop over all neighbors (this now includes the cell itself!) to iterate
  // over all beads of type2 of the cell and its neighbors
  for (vector<cell_t *>::iterator iterc2 = cell._neighbours.begin();
       iterc2 != cell._neighbours.end(); ++iterc2) {
    for (iter2 = (*(*iterc2))._beads2.begin();
         iter2 != (*(*iterc2))._beads2.end(); ++iter2) {

      if (bead == *iter2) {
        continue;
      }

      // loop again over all neighbors (this now includes the cell itself!)
      // to iterate over all beads of type3 of the cell and its neighbors
      for (auto &_neighbour : cell._neighbours) {
        for (iter3 = (*_neighbour)._beads3.begin();
             iter3 != (*_neighbour)._beads3.end(); ++iter3) {

          // do not include the same beads twice in one triple!
          if (bead == *iter3) {
            continue;
          }
          if (*iter2 == *iter3) {
            continue;
          }

          Eigen::Vector3d v = (*iter2)->getPos();
          Eigen::Vector3d z = (*iter3)->getPos();

          Eigen::Vector3d r12 = _top->BCShortestConnection(u, v);
          Eigen::Vector3d r13 = _top->BCShortestConnection(u, z);
          Eigen::Vector3d r23 = _top->BCShortestConnection(v, z);
          double d12 = r12.norm();
          double d13 = r13.norm();
          double d23 = r23.norm();

          // to do: at the moment use only one cutoff value
          // to do: so far only check the distance between bead 1 (central
          // bead) and bead2 and bead 3
          if ((d12 < _cutoff) && (d13 < _cutoff)) {
            /// experimental: at the moment exclude interaction as soon as
            /// one of the three pairs (1,2) (1,3) (2,3) is excluded!
            if (_do_exclusions) {
              if ((_top->getExclusions().IsExcluded(bead, *iter2)) ||
                  (_top->getExclusions().IsExcluded(bead, *iter3)) ||
                  (_top->getExclusions().IsExcluded(*iter2, *iter3))) {
                continue;
              }
            }
            if ((*_match_function)(bead, *iter2, *iter3, r12, r13, r23, d12,
                                   d13, d23)) {
              if (!FindTriple(bead, *iter2, *iter3)) {
                AddTriple(_triple_creator(bead, *iter2, *iter3, r12, r13, r23));
              }
            }
          }
        }
      }
    }
  }
}

}  // namespace csg
}  // namespace votca
