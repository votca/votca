/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_REGULAR_GRID_H
#define VOTCA_XTP_REGULAR_GRID_H

// Local VOTCA includes
#include "gridbox.h"

namespace votca {
namespace xtp {
class QMMolecule;
class aobasis;

class Regular_Grid {
 public:
  void GridSetup(const Eigen::Array<Index, 3, 1>& steps,
                 const Eigen::Array3d& padding, const QMMolecule& atoms,
                 const AOBasis& basis);

  void GridSetup(const Eigen::Array3d& stepsize, const Eigen::Array3d& padding,
                 const QMMolecule& atoms, const AOBasis& basis);

  Index getGridSize() const { return totalgridsize_; }
  Index getBoxesSize() const { return Index(grid_boxes_.size()); }

  const GridBox& operator[](Index index) const { return grid_boxes_[index]; }
  GridBox& operator[](Index index) { return grid_boxes_[index]; }

  std::vector<GridBox>::iterator begin() { return grid_boxes_.begin(); }
  std::vector<GridBox>::iterator end() { return grid_boxes_.end(); }

  std::vector<GridBox>::const_iterator begin() const {
    return grid_boxes_.begin();
  }
  std::vector<GridBox>::const_iterator end() const { return grid_boxes_.end(); }

  Eigen::Array3d getStepSizes() const { return stepsizes_; }

  Eigen::Vector3d getStartingPoint() const { return startingpoint_; }

  Eigen::Array<Index, 3, 1> getSteps() const { return steps_; }

 private:
  Index totalgridsize_;
  std::vector<GridBox> grid_boxes_;

  Eigen::Array3d stepsizes_ = Eigen::Array3d::Zero();
  Eigen::Vector3d startingpoint_ = Eigen::Vector3d::Zero();
  Eigen::Array<Index, 3, 1> steps_ = Eigen::Array<Index, 3, 1>::Zero();
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_REGULAR_GRID_H
