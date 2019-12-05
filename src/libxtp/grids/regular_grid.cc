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

#include <votca/xtp/qmmolecule.h>
#include <votca/xtp/radial_euler_maclaurin_rule.h>
#include <votca/xtp/regular_grid.h>
#include <votca/xtp/sphere_lebedev_rule.h>

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
  Eigen::Array3d minpos = min - padding;
  Eigen::Array3d stepsizes = steps.cast<double>() / (max - min + 2 * padding);

  const Index boxsize = 10;

  Eigen::Array<Index, 3, 1> tempgridsize =
      (steps.cast<double>() / boxsize).ceil().cast<Index>();

  std::vector<std::vector<
      std::vector<std::vector<GridContainers::Cartesian_gridpoint>>>>
      tempgrid = std::vector<std::vector<
          std::vector<std::vector<GridContainers::Cartesian_gridpoint>>>>(
          tempgridsize.x());
  for (Index i = 0; i < tempgridsize.x(); i++) {
    tempgrid[i].resize(tempgridsize.y());
    for (Index j = 0; j < tempgridsize.y(); j++) {
      tempgrid[i][j].resize(tempgridsize.z());
      for (Index k = 0; k < tempgridsize.z(); k++) {
        tempgrid[i][j][k].reserve(std::pow(boxsize, 3));
      }
    }
  }

  for (Index i = 0; i < steps.x(); i++) {
    double x = minpos.x() + double(i) * stepsizes.x();
    Index box_index_x = i / tempgridsize.x();
    for (Index j = 0; j < steps.y(); j++) {
      double y = minpos.y() + double(j) * stepsizes.y();
      Index box_index_y = j / tempgridsize.y();
      for (Index k = 0; k < steps.z(); k++) {
        double z = minpos.z() + double(k) * stepsizes.z();
        Index box_index_z = k / tempgridsize.z();
        GridContainers::Cartesian_gridpoint point;
        point.grid_weight = 1.0;
        point.grid_pos = Eigen::Vector3d(x, y, z);
        tempgrid[box_index_x][box_index_y][box_index_z].push_back(point);
      }
    }
  }

  for (auto& boxes_xy : tempgrid) {
    for (auto& boxes_z : boxes_xy) {
      for (auto& box : boxes_z) {
        if (box.empty()) {
          continue;
        }
        GridBox gridbox;
        for (const auto& point : box) {
          gridbox.addGridPoint(point);
        }
        _grid_boxes.push_back(gridbox);
      }
    }
  }

#pragma omp parallel for
  for (Index i = 0; i < getBoxesSize(); i++) {
    _grid_boxes[i].FindSignificantShells(basis);
  }

  _totalgridsize = 0;
  for (auto& box : _grid_boxes) {
    _totalgridsize += box.size();
    box.PrepareForIntegration();
  }
}

}  // namespace xtp
}  // namespace votca
