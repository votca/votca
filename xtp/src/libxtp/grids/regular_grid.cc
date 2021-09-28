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

// Local VOTCA includes
#include "votca/xtp/regular_grid.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/radial_euler_maclaurin_rule.h"
#include "votca/xtp/sphere_lebedev_rule.h"

namespace votca {
namespace xtp {

void Regular_Grid::GridSetup(const Eigen::Array3d& stepsizes,
                             const Eigen::Array3d& padding,
                             const QMMolecule& atoms, const AOBasis& basis) {

  std::pair<Eigen::Vector3d, Eigen::Vector3d> extension =
      atoms.CalcSpatialMinMax();
  Eigen::Array3d min = extension.first.array();
  Eigen::Array3d max = extension.second.array();
  Eigen::Array3d doublesteps = (max - min + 2 * padding) / stepsizes + 1.0;
  Eigen::Array<votca::Index, 3, 1> steps = (doublesteps.ceil()).cast<Index>();

  // needed to symmetrize grid around molecule
  Eigen::Array3d padding_sym =
      (doublesteps - steps.cast<double>()) * stepsizes * 0.5 + padding;
  GridSetup(steps, padding_sym, atoms, basis);
}

void Regular_Grid::GridSetup(const Eigen::Array<Index, 3, 1>& steps,
                             const Eigen::Array3d& padding,
                             const QMMolecule& atoms, const AOBasis& basis) {

  std::pair<Eigen::Vector3d, Eigen::Vector3d> extension =
      atoms.CalcSpatialMinMax();
  Eigen::Array3d min = extension.first.array();
  Eigen::Array3d max = extension.second.array();
  startingpoint_ = min - padding;
  stepsizes_ = (max - min + 2 * padding) / (steps - 1).cast<double>();
  steps_ = steps;
  const Index gridboxsize = 500;

  GridBox gridbox;
  for (Index i = 0; i < steps.x(); i++) {
    double x = startingpoint_.x() + double(i) * stepsizes_.x();
    for (Index j = 0; j < steps.y(); j++) {
      double y = startingpoint_.y() + double(j) * stepsizes_.y();
      for (Index k = 0; k < steps.z(); k++) {
        double z = startingpoint_.z() + double(k) * stepsizes_.z();
        GridContainers::Cartesian_gridpoint point;
        point.grid_weight = 1.0;
        point.grid_pos = Eigen::Vector3d(x, y, z);
        gridbox.addGridPoint(point);
        if (gridbox.size() == gridboxsize) {
          grid_boxes_.push_back(gridbox);
          gridbox = GridBox();
        }
      }
    }
  }
  if (gridbox.size() > 0) {
    grid_boxes_.push_back(gridbox);
  }
#pragma omp parallel for
  for (Index i = 0; i < getBoxesSize(); i++) {
    grid_boxes_[i].FindSignificantShells(basis);
    grid_boxes_[i].PrepareForIntegration();
  }

  totalgridsize_ = 0;
  for (auto& box : grid_boxes_) {
    totalgridsize_ += box.size();
  }
}

}  // namespace xtp
}  // namespace votca
