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
#include "votca/csg/nblistgrid.h"
#include "votca/csg/topology.h"
#include "votca/tools/NDimVector.h"

namespace votca {
namespace csg {

using namespace std;

void NBListGrid::Generate(BeadList &list1, BeadList &list2,
                          bool do_exclusions) {

  do_exclusions_ = do_exclusions;
  if (list1.empty()) {
    return;
  }
  if (list2.empty()) {
    return;
  }

  assert(&(list1.getTopology()) == &(list2.getTopology()));
  const Topology &top = list1.getTopology();

  InitializeGrid(top.getBox());

  // Add all beads of list1
  for (auto &iter : list1) {
    getCell(iter->getPos()).beads_.push_back(iter);
  }

  for (auto &iter : list2) {
    cell_t &cell = getCell(iter->getPos());
    TestBead(top, cell, iter);
  }
}

void NBListGrid::Generate(BeadList &list, bool do_exclusions) {
  do_exclusions_ = do_exclusions;
  if (list.empty()) {
    return;
  }

  const Topology &top = list.getTopology();

  InitializeGrid(top.getBox());

  for (auto &iter : list) {
    cell_t &cell = getCell(iter->getPos());
    TestBead(top, cell, iter);
    getCell(iter->getPos()).beads_.push_back(iter);
  }
}
void NBListGrid::InitializeGrid(const Eigen::Matrix3d &box) {
  box_a_ = box.col(0);
  box_b_ = box.col(1);
  box_c_ = box.col(2);

  // create plane normals
  norm_a_ = box_b_.cross(box_c_);
  norm_b_ = box_c_.cross(box_a_);
  norm_c_ = box_a_.cross(box_b_);

  norm_a_.normalize();
  norm_b_.normalize();
  norm_c_.normalize();

  double la = box_a_.dot(norm_a_);
  double lb = box_b_.dot(norm_b_);
  double lc = box_c_.dot(norm_c_);

  // calculate grid size, each grid has to be at least size of cut-off
  box_Na_ = Index(std::max(std::abs(la / cutoff_), 1.0));
  box_Nb_ = Index(std::max(std::abs(lb / cutoff_), 1.0));
  box_Nc_ = Index(std::max(std::abs(lc / cutoff_), 1.0));

  norm_a_ = norm_a_ / box_a_.dot(norm_a_) * (double)box_Na_;
  norm_b_ = norm_b_ / box_b_.dot(norm_b_) * (double)box_Nb_;
  norm_c_ = norm_c_ / box_c_.dot(norm_c_) * (double)box_Nc_;

  grid_ = tools::NDimVector<cell_t, 3>(box_Na_, box_Nb_, box_Nc_);

  Index a1, a2, b1, b2, c1, c2;

  a1 = b1 = c1 = -1;
  a2 = b2 = c2 = 1;

  if (box_Na_ < 3) {
    a2 = 0;
  }
  if (box_Nb_ < 3) {
    b2 = 0;
  }
  if (box_Nc_ < 3) {
    c2 = 0;
  }

  if (box_Na_ < 2) {
    a1 = 0;
  }
  if (box_Nb_ < 2) {
    b1 = 0;
  }
  if (box_Nc_ < 2) {
    c1 = 0;
  }

  // wow, setting up the neighbours is an ugly for construct!
  // loop from N..2*N to avoid if and only use %
  for (Index a = box_Na_; a < 2 * box_Na_; ++a) {
    for (Index b = box_Nb_; b < 2 * box_Nb_; ++b) {
      for (Index c = box_Nc_; c < 2 * box_Nc_; ++c) {
        cell_t &cell = grid_(a % box_Na_, b % box_Nb_, c % box_Nc_);
        for (Index aa = a + a1; aa <= a + a2; ++aa) {
          for (Index bb = b + b1; bb <= b + b2; ++bb) {
            for (Index cc = c + c1; cc <= c + c2; ++cc) {
              cell_t *cell2 = &grid_(aa % box_Na_, bb % box_Nb_, cc % box_Nc_);
              if (cell2 == &cell) {
                continue;  // ignore self
              }
              cell.neighbours_.push_back(
                  &grid_(aa % box_Na_, bb % box_Nb_, cc % box_Nc_));
            }
          }
        }
      }
    }
  }
}

NBListGrid::cell_t &NBListGrid::getCell(const Eigen::Vector3d &r) {
  Index a = (Index)floor(r.dot(norm_a_));
  Index b = (Index)floor(r.dot(norm_b_));
  Index c = (Index)floor(r.dot(norm_c_));

  if (a < 0) {
    a = box_Na_ + a % box_Na_;
  }
  a %= box_Na_;

  if (b < 0) {
    b = box_Nb_ + b % box_Nb_;
  }
  b %= box_Nb_;

  if (c < 0) {
    c = box_Nc_ + c % box_Nc_;
  }
  c %= box_Nc_;

  return grid_(a, b, c);
}

void NBListGrid::TestBead(const Topology &top, NBListGrid::cell_t &cell,
                          Bead *bead) {
  TestCell(top, cell, bead);
  for (auto &neighbour : cell.neighbours_) {
    TestCell(top, *neighbour, bead);
  }
}

void NBListGrid::TestCell(const Topology &top, NBListGrid::cell_t &cell,
                          Bead *bead) {
  const Eigen::Vector3d &u = bead->getPos();

  for (auto &bead_ : cell.beads_) {

    const Eigen::Vector3d &v = bead_->getPos();
    const Eigen::Vector3d &r = top.BCShortestConnection(v, u);
    double d = r.norm();
    if (d < cutoff_) {
      if (do_exclusions_) {
        if (top.getExclusions().IsExcluded(bead_, bead)) {
          continue;
        }
      }
      if ((*match_function_)(bead_, bead, r, d)) {
        if (!FindPair(bead_, bead)) {
          AddPair(pair_creator_(bead_, bead, r));
        }
      }
    }
  }
}

}  // namespace csg
}  // namespace votca
