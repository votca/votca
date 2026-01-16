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
#ifndef VOTCA_XTP_GRIDBOX_H
#define VOTCA_XTP_GRIDBOX_H

// Local VOTCA includes
#include "aoshell.h"
#include "grid_containers.h"

namespace votca {
namespace xtp {

struct GridboxRange {
  Index start;
  Index size;
};
class GridBox {

 public:
  void FindSignificantShells(const AOBasis& basis);
  AOShell::AOValues CalcAOValues(const Eigen::Vector3d& point) const;

  const std::vector<Eigen::Vector3d>& getGridPoints() const { return grid_pos; }

  const std::vector<double>& getGridWeights() const { return weights; }
  std::vector<double>& getPotentialValues()  { return potential_values; }

  const std::vector<const AOShell*>& getShells() const {
    return significant_shells;
  }

  const std::vector<GridboxRange>& getAOranges() const { return aoranges; }

  Index size() const { return Index(grid_pos.size()); }

  Index Shellsize() const { return Index(significant_shells.size()); }

  Index Matrixsize() const { return matrix_size; }

  void addGridBox(const GridBox& box) {
    grid_pos.insert(grid_pos.end(), box.grid_pos.begin(), box.grid_pos.end());
    weights.insert(weights.end(), box.weights.begin(), box.weights.end());
    potential_values.insert(potential_values.end(), box.potential_values.begin(), box.potential_values.end());
    return;
  }

  void addGridPoint(const GridContainers::Cartesian_gridpoint& point) {
    grid_pos.push_back(point.grid_pos);
    weights.push_back(point.grid_weight);
    potential_values.push_back(0.0);
  };

  void addShell(const AOShell* shell) {
    significant_shells.push_back(shell);
    matrix_size += shell->getNumFunc();
  };

  void PrepareForIntegration();

  Eigen::MatrixXd ReadFromBigMatrix(const Eigen::MatrixXd& bigmatrix) const;

  Eigen::VectorXd ReadFromBigVector(const Eigen::VectorXd& bigvector) const;

  void AddtoBigMatrix(Eigen::MatrixXd& bigmatrix,
                      const Eigen::MatrixXd& smallmatrix) const;

  static bool compareGridboxes(GridBox& box1, GridBox& box2) {
    if (box1.Matrixsize() != box2.Matrixsize()) {
      return false;
    }
    if (box1.Shellsize() != box2.Shellsize()) {
      return false;
    }
    for (Index i = 0; i < Index(box1.significant_shells.size()); ++i) {
      if (box1.significant_shells[i] != box2.significant_shells[i]) {
        return false;
      }
    }
    return true;
  }

 private:
  Index matrix_size = 0;
  std::vector<GridboxRange> aoranges;
  std::vector<GridboxRange> ranges;
  std::vector<GridboxRange> inv_ranges;
  std::vector<Eigen::Vector3d> grid_pos;
  std::vector<const AOShell*> significant_shells;
  std::vector<double> weights;
  std::vector<double> potential_values;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_GRIDBOX_H
