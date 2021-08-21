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
#include "votca/csg/nblistgrid_3body.h"
#include "votca/csg/topology.h"

namespace votca {
namespace csg {

using namespace std;

void NBListGrid_3Body::Generate(BeadList &list1, BeadList &list2,
                                BeadList &list3, bool do_exclusions) {

  do_exclusions_ = do_exclusions;
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
  const Topology &top = list1.getTopology();

  InitializeGrid(top.getBox());

  // Add all beads of list1 to  beads1_
  for (auto &iter : list1) {
    getCell(iter->getPos()).beads1_.push_back(iter);
  }

  // Add all beads of list2 to  beads2_
  for (auto &iter : list2) {
    getCell(iter->getPos()).beads2_.push_back(iter);
  }

  // Add all beads of list2 to  beads3_
  for (auto &iter : list3) {
    getCell(iter->getPos()).beads3_.push_back(iter);
  }

  // loop over beads of list 1 again to get the correlations
  for (auto &iter : list1) {
    cell_t &cell = getCell(iter->getPos());
    TestBead(top, cell, iter);
  }
}

void NBListGrid_3Body::Generate(BeadList &list1, BeadList &list2,
                                bool do_exclusions) {
  do_exclusions_ = do_exclusions;
  if (list1.empty()) {
    return;
  }
  if (list2.empty()) {
    return;
  }

  // check if both bead lists "have" the same topology
  assert(&(list1.getTopology()) == &(list2.getTopology()));
  const Topology &top = list1.getTopology();

  InitializeGrid(top.getBox());

  // Add all beads of list1 to  beads1_
  for (auto &bead : list1) {
    getCell(bead->getPos()).beads1_.push_back(bead);
  }

  // Add all beads of list2 to  beads2_
  for (auto &bead : list2) {
    getCell(bead->getPos()).beads2_.push_back(bead);
  }

  // In this case type2 and type3 are the same
  for (auto &cell : grid_) {
    cell.beads3_ = cell.beads2_;
  }

  // loop over beads of list 1 again to get the correlations
  for (auto &bead : list1) {
    cell_t &cell = getCell(bead->getPos());
    TestBead(top, cell, bead);
  }
}

void NBListGrid_3Body::Generate(BeadList &list, bool do_exclusions) {
  do_exclusions_ = do_exclusions;
  if (list.empty()) {
    return;
  }

  const Topology &top = list.getTopology();

  InitializeGrid(top.getBox());

  // Add all beads of list to all! bead lists of the cell
  for (auto &iter : list) {
    getCell(iter->getPos()).beads1_.push_back(iter);
  }

  for (auto &cell : grid_) {
    cell.beads2_ = cell.beads1_;
    cell.beads3_ = cell.beads1_;
  }

  // loop over beads again to get the correlations (as all of the same type
  // here)
  for (auto &bead : list) {
    cell_t &cell = getCell(bead->getPos());
    TestBead(top, cell, bead);
  }
}

void NBListGrid_3Body::InitializeGrid(const Eigen::Matrix3d &box) {
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

  norm_a_ = norm_a_ / la * (double)box_Na_;
  norm_b_ = norm_b_ / lb * (double)box_Nb_;
  norm_c_ = norm_c_ / lc * (double)box_Nc_;

  grid_.resize(box_Na_ * box_Nb_ * box_Nc_);

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
        cell_t &cell = getCell(a % box_Na_, b % box_Nb_, c % box_Nc_);
        for (Index aa = a + a1; aa <= a + a2; ++aa) {
          for (Index bb = b + b1; bb <= b + b2; ++bb) {
            for (Index cc = c + c1; cc <= c + c2; ++cc) {
              // test: for 3body algorithm: each cell is a neighbor of its own
              // !!!
              cell.neighbours_.push_back(
                  &getCell(aa % box_Na_, bb % box_Nb_, cc % box_Nc_));
            }
          }
        }
      }
    }
  }
}

NBListGrid_3Body::cell_t &NBListGrid_3Body::getCell(const Eigen::Vector3d &r) {
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

  return getCell(a, b, c);
}

void NBListGrid_3Body::TestBead(const Topology &top,
                                NBListGrid_3Body::cell_t &cell, Bead *bead) {
  BeadList::iterator iter2;
  BeadList::iterator iter3;
  Eigen::Vector3d u = bead->getPos();

  // loop over all neighbors (this now includes the cell itself!) to iterate
  // over all beads of type2 of the cell and its neighbors
  for (vector<cell_t *>::iterator iterc2 = cell.neighbours_.begin();
       iterc2 != cell.neighbours_.end(); ++iterc2) {
    for (iter2 = (*(*iterc2)).beads2_.begin();
         iter2 != (*(*iterc2)).beads2_.end(); ++iter2) {

      if (bead == *iter2) {
        continue;
      }

      // loop again over all neighbors (this now includes the cell itself!)
      // to iterate over all beads of type3 of the cell and its neighbors
      for (auto &neighbour_ : cell.neighbours_) {
        for (iter3 = (*neighbour_).beads3_.begin();
             iter3 != (*neighbour_).beads3_.end(); ++iter3) {

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
          if ((d12 < cutoff_) && (d13 < cutoff_)) {
            /// experimental: at the moment exclude interaction as soon as
            /// one of the three pairs (1,2) (1,3) (2,3) is excluded!
            if (do_exclusions_) {
              if ((top.getExclusions().IsExcluded(bead, *iter2)) ||
                  (top.getExclusions().IsExcluded(bead, *iter3)) ||
                  (top.getExclusions().IsExcluded(*iter2, *iter3))) {
                continue;
              }
            }
            if ((*match_function_)(bead, *iter2, *iter3, r12, r13, r23, d12,
                                   d13, d23)) {
              if (!FindTriple(bead, *iter2, *iter3)) {
                AddTriple(triple_creator_(bead, *iter2, *iter3, r12, r13, r23));
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
