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
  assert(&(list1.getTopology()) == &(list2.getTopology()));
  assert(&(list1.getTopology()) == &(list3.getTopology()));
  assert(&(list2.getTopology()) == &(list3.getTopology()));
  const Topology & top = list1.getTopology();

  InitializeGrid(top.getBox());

  // Add all beads of list1 to _beads1
  for (auto &iter : list1) {
    getCell(iter->getPos())._beads1.push_back(iter);
  }

  // Add all beads of list2 to _beads2
  for (auto &iter : list2) {
    getCell(iter->getPos())._beads2.push_back(iter);
  }

  // Add all beads of list2 to _beads3
  for (auto &iter : list3) {
    getCell(iter->getPos())._beads3.push_back(iter);
  }

  // loop over beads of list 1 again to get the correlations
  for (auto &iter : list1) {
    cell_t &cell = getCell(iter->getPos());
    TestBead(top, cell, iter);
  }
}

void NBListGrid_3Body::Generate(BeadList &list1, BeadList &list2,
                                bool do_exclusions) {
  _do_exclusions = do_exclusions;
  if (list1.empty()) {
    return;
  }
  if (list2.empty()) {
    return;
  }

  // check if both bead lists "have" the same topology
  assert(list1.getTopology() == list2.getTopology());
  const Topology & top = list1.getTopology();

  InitializeGrid(top.getBox());

  // Add all beads of list1 to _beads1
  for (auto &bead : list1) {
    getCell(bead->getPos())._beads1.push_back(bead);
  }

  // Add all beads of list2 to _beads2
  for (auto &bead : list2) {
    getCell(bead->getPos())._beads2.push_back(bead);
  }

  // In this case type2 and type3 are the same
  for (auto &cell : _grid) {
    cell._beads3 = cell._beads2;
  }

  // loop over beads of list 1 again to get the correlations
  for (auto &bead : list1) {
    cell_t &cell = getCell(bead->getPos());
    TestBead(top, cell, bead);
  }
}

void NBListGrid_3Body::Generate(BeadList &list, bool do_exclusions) {
  _do_exclusions = do_exclusions;
  if (list.empty()) {
    return;
  }

  const Topology & top = list.getTopology();

  InitializeGrid(top.getBox());

  // Add all beads of list to all! bead lists of the cell
  for (auto &iter : list) {
    getCell(iter->getPos())._beads1.push_back(iter);
  }

  for (auto &cell : _grid) {
    cell._beads2 = cell._beads1;
    cell._beads3 = cell._beads1;
  }

  // loop over beads again to get the correlations (as all of the same type
  // here)
  for (auto &bead : list) {
    cell_t &cell = getCell(bead->getPos());
    TestBead(top, cell, bead);
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
  _box_Na = Index(std::max(std::abs(la / _cutoff), 1.0));
  _box_Nb = Index(std::max(std::abs(lb / _cutoff), 1.0));
  _box_Nc = Index(std::max(std::abs(lc / _cutoff), 1.0));

  _norm_a = _norm_a / la * (double)_box_Na;
  _norm_b = _norm_b / lb * (double)_box_Nb;
  _norm_c = _norm_c / lc * (double)_box_Nc;

  _grid.resize(_box_Na * _box_Nb * _box_Nc);

  Index a1, a2, b1, b2, c1, c2;

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
  for (Index a = _box_Na; a < 2 * _box_Na; ++a) {
    for (Index b = _box_Nb; b < 2 * _box_Nb; ++b) {
      for (Index c = _box_Nc; c < 2 * _box_Nc; ++c) {
        cell_t &cell = getCell(a % _box_Na, b % _box_Nb, c % _box_Nc);
        for (Index aa = a + a1; aa <= a + a2; ++aa) {
          for (Index bb = b + b1; bb <= b + b2; ++bb) {
            for (Index cc = c + c1; cc <= c + c2; ++cc) {
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
  Index a = (Index)floor(r.dot(_norm_a));
  Index b = (Index)floor(r.dot(_norm_b));
  Index c = (Index)floor(r.dot(_norm_c));

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

void NBListGrid_3Body::TestBead(const Topology & top, NBListGrid_3Body::cell_t &cell, Bead *bead) {
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

          Eigen::Vector3d r12 = top.BCShortestConnection(u, v);
          Eigen::Vector3d r13 = top.BCShortestConnection(u, z);
          Eigen::Vector3d r23 = top.BCShortestConnection(v, z);
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
              if ((top.getExclusions().IsExcluded(bead, *iter2)) ||
                  (top.getExclusions().IsExcluded(bead, *iter3)) ||
                  (top.getExclusions().IsExcluded(*iter2, *iter3))) {
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
