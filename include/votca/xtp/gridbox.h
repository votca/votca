/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef __XTP_GRIDBOX__H
#define __XTP_GRIDBOX__H

#include <votca/tools/vec.h>
#include <votca/xtp/aoshell.h>
#include <votca/xtp/grid_containers.h>

namespace votca {
namespace xtp {

struct GridboxRange {
  int start;
  int size;
};
class GridBox {

 public:
  const std::vector<tools::vec>& getGridPoints() const { return grid_pos; }

  const std::vector<double>& getGridWeights() const { return weights; }

  const std::vector<const AOShell*>& getShells() const {
    return significant_shells;
  }

  const std::vector<GridboxRange>& getAOranges() const { return aoranges; }

  unsigned size() const { return grid_pos.size(); }

  unsigned Shellsize() const { return significant_shells.size(); }

  unsigned Matrixsize() const { return matrix_size; }

  void addGridBox(const GridBox& box) {
    const std::vector<tools::vec>& p = box.getGridPoints();
    const std::vector<double>& w = box.getGridWeights();
    for (unsigned i = 0; i < w.size(); ++i) {
      grid_pos.push_back(p[i]);
      weights.push_back(w[i]);
    }
    return;
  }

  void addGridPoint(const GridContainers::Cartesian_gridpoint& point) {
    grid_pos.push_back(point.grid_pos);
    weights.push_back(point.grid_weight);
  };

  void addShell(const AOShell* shell) {
    significant_shells.push_back(shell);
    matrix_size += shell->getNumFunc();
  };

  void prepareDensity() { densities.reserve(grid_pos.size()); }

  void addDensity(double density) { densities.push_back(density); }

  const std::vector<double>& getGridDensities() const { return densities; }

  void PrepareForIntegration();

  Eigen::MatrixXd ReadFromBigMatrix(const Eigen::MatrixXd& bigmatrix) const;

  void AddtoBigMatrix(Eigen::MatrixXd& bigmatrix,
                      const Eigen::MatrixXd& smallmatrix) const;

  void setIndexoffirstgridpoint(unsigned indexoffirstgridpoint) {
    _indexoffirstgridpoint = indexoffirstgridpoint;
  }
  unsigned getIndexoffirstgridpoint() const { return _indexoffirstgridpoint; }

  static bool compareGridboxes(GridBox& box1, GridBox& box2) {
    if (box1.Matrixsize() != box2.Matrixsize()) {
      return false;
    }
    if (box1.Shellsize() != box2.Shellsize()) {
      return false;
    }
    for (unsigned i = 0; i < box1.significant_shells.size(); ++i) {
      if (box1.significant_shells[i] != box2.significant_shells[i]) {
        return false;
      }
    }
    return true;
  }

 private:
  unsigned _indexoffirstgridpoint;
  unsigned matrix_size = 0;
  std::vector<GridboxRange> aoranges;
  std::vector<GridboxRange> ranges;
  std::vector<GridboxRange> inv_ranges;
  std::vector<tools::vec> grid_pos;  // bohr
  std::vector<const AOShell*> significant_shells;
  std::vector<double> weights;
  std::vector<double> densities;
  std::vector<Eigen::MatrixXd> dens_grad;
};

}  // namespace xtp
}  // namespace votca
#endif /* NUMERICAL_INTEGRATION_H */
