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

#include <votca/csg/nblistgrid.h>

namespace votca {
namespace csg {

using namespace std;

void NBListGrid::Generate(BeadList &list1, BeadList &list2,
                          bool do_exclusions) {
  BeadList::iterator iter;
  _do_exclusions = do_exclusions;
  if (list1.empty()) return;
  if (list2.empty()) return;

  assert(list1.getTopology() == list2.getTopology());
  Topology *top = _top = list1.getTopology();

  InitializeGrid(top->getBox());

  // Add all beads of list1
  for (iter = list1.begin(); iter != list1.end(); ++iter)
    getCell((*iter)->getPos())._beads.push_back(*iter);

  for (iter = list2.begin(); iter != list2.end(); ++iter) {
    cell_t &cell = getCell((*iter)->getPos());
    TestBead(cell, *iter);
  }
}

void NBListGrid::Generate(BeadList &list, bool do_exclusions) {
  BeadList::iterator iter;
  _do_exclusions = do_exclusions;
  if (list.empty()) return;

  Topology *top = _top = list.getTopology();

  InitializeGrid(top->getBox());

  for (iter = list.begin(); iter != list.end(); ++iter) {
    cell_t &cell = getCell((*iter)->getPos());
    TestBead(cell, *iter);
    getCell((*iter)->getPos())._beads.push_back(*iter);
  }
}
void NBListGrid::InitializeGrid(const Eigen::Matrix3d &box) {
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

  _norm_a = _norm_a / _box_a.dot(_norm_a) * (double)_box_Na;
  _norm_b = _norm_b / _box_b.dot(_norm_b) * (double)_box_Nb;
  _norm_c = _norm_c / _box_c.dot(_norm_c) * (double)_box_Nc;

  _grid.resize(_box_Na * _box_Nb * _box_Nc);

  int a1, a2, b1, b2, c1, c2;

  a1 = b1 = c1 = -1;
  a2 = b2 = c2 = 1;

  if (_box_Na < 3) a2 = 0;
  if (_box_Nb < 3) b2 = 0;
  if (_box_Nc < 3) c2 = 0;

  if (_box_Na < 2) a1 = 0;
  if (_box_Nb < 2) b1 = 0;
  if (_box_Nc < 2) c1 = 0;

  // wow, setting up the neighbours is an ugly for construct!
  // loop from N..2*N to avoid if and only use %
  for (int a = _box_Na; a < 2 * _box_Na; ++a)
    for (int b = _box_Nb; b < 2 * _box_Nb; ++b)
      for (int c = _box_Nc; c < 2 * _box_Nc; ++c) {
        cell_t &cell = getCell(a % _box_Na, b % _box_Nb, c % _box_Nc);
        for (int aa = a + a1; aa <= a + a2; ++aa)
          for (int bb = b + b1; bb <= b + b2; ++bb)
            for (int cc = c + c1; cc <= c + c2; ++cc) {
              cell_t *c = &getCell(aa % _box_Na, bb % _box_Nb, cc % _box_Nc);
              if (c == &cell) continue;  // ignore self
              cell._neighbours.push_back(
                  &getCell(aa % _box_Na, bb % _box_Nb, cc % _box_Nc));
            }
      }
}

NBListGrid::cell_t &NBListGrid::getCell(const Eigen::Vector3d &r) {
  int a = (int)floor(r.dot(_norm_a));
  int b = (int)floor(r.dot(_norm_b));
  int c = (int)floor(r.dot(_norm_c));

  if (a < 0) a = _box_Na + a % _box_Na;
  a %= _box_Na;

  if (b < 0) b = _box_Nb + b % _box_Nb;
  b %= _box_Nb;

  if (c < 0) c = _box_Nc + c % _box_Nc;
  c %= _box_Nc;

  return getCell(a, b, c);
}

void NBListGrid::TestBead(NBListGrid::cell_t &cell, Bead *bead) {
  TestCell(cell, bead);
  for (vector<cell_t *>::iterator iter = cell._neighbours.begin();
       iter != cell._neighbours.end(); ++iter) {
    TestCell(*(*iter), bead);
  }
}

void NBListGrid::TestCell(NBListGrid::cell_t &cell, Bead *bead) {
  BeadList::iterator iter;
  Eigen::Vector3d u = bead->getPos();

  for (iter = cell._beads.begin(); iter != cell._beads.end(); ++iter) {

    Eigen::Vector3d v = (*iter)->getPos();
    Eigen::Vector3d r = _top->BCShortestConnection(v, u);
    double d = r.norm();
    if (d < _cutoff) {
      if (_do_exclusions)
        if (_top->getExclusions().IsExcluded((*iter), bead)) {
          continue;
        }
      if ((*_match_function)(*iter, bead, r, d))
        if (!FindPair(*iter, bead)) AddPair(_pair_creator(*iter, bead, r));
    }
  }
}

}  // namespace csg
}  // namespace votca
