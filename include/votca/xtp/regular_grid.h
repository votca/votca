/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#ifndef XTP_REGULAR_GRID_H
#define XTP_REGULAR_GRID_H

#include <votca/xtp/gridbox.h>

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

  Index getGridSize() const { return _totalgridsize; }
  Index getBoxesSize() const { return Index(_grid_boxes.size()); }

  const GridBox& operator[](Index index) const { return _grid_boxes[index]; }
  GridBox& operator[](Index index) { return _grid_boxes[index]; }

  std::vector<GridBox>::iterator begin() { return _grid_boxes.begin(); }
  std::vector<GridBox>::iterator end() { return _grid_boxes.end(); }

  std::vector<GridBox>::const_iterator begin() const {
    return _grid_boxes.begin();
  }
  std::vector<GridBox>::const_iterator end() const { return _grid_boxes.end(); }

  Eigen::Array3d getStepSizes() const { return _stepsizes; }

  Eigen::Vector3d getStartingPoint() const { return _startingpoint; }

  Eigen::Array<Index, 3, 1> getSteps() const { return _steps; }

 private:
  Index _totalgridsize;
  std::vector<GridBox> _grid_boxes;

  Eigen::Array3d _stepsizes = Eigen::Array3d::Zero();
  Eigen::Vector3d _startingpoint = Eigen::Vector3d::Zero();
  Eigen::Array<Index, 3, 1> _steps = Eigen::Array<Index, 3, 1>::Zero();
};

}  // namespace xtp
}  // namespace votca
#endif  // XTP_VXC_GRID_H
