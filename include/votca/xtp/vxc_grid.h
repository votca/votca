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
#ifndef XTP_VXC_GRID_H
#define XTP_VXC_GRID_H

#include <votca/xtp/grid_containers.h>
#include <votca/xtp/gridbox.h>

namespace votca {
namespace xtp {
class LebedevGrid;
class QMMolecule;
class aobasis;

class Vxc_Grid {
 public:
  void GridSetup(const std::string& type, const QMMolecule& atoms,
                 const AOBasis& basis);

  std::vector<const Eigen::Vector3d*> getGridpoints() const;
  std::vector<double> getWeightedDensities() const;
  Index getGridSize() const { return _totalgridsize; }
  Index getBoxesSize() const { return Index(_grid_boxes.size()); }

  const GridBox& operator[](Index index) const { return _grid_boxes[index]; }
  GridBox& operator[](Index index) { return _grid_boxes[index]; }

  typename std::vector<GridBox>::iterator begin() {
    return _grid_boxes.begin();
  }
  typename std::vector<GridBox>::iterator end() { return _grid_boxes.end(); }

  typename std::vector<GridBox>::const_iterator begin() const {
    return _grid_boxes.begin();
  }
  typename std::vector<GridBox>::const_iterator end() const {
    return _grid_boxes.end();
  }

 private:
  void FindSignificantShells(const AOBasis& basis);

  double erf1c(double x) const;

  void SortGridpointsintoBlocks(
      const std::vector<std::vector<GridContainers::Cartesian_gridpoint> >&
          grid);

  Eigen::MatrixXd CalcInverseAtomDist(const QMMolecule& atoms) const;
  Index UpdateOrder(LebedevGrid& sphericalgridofElement, Index maxorder,
                    std::vector<double>& PruningIntervals, double r) const;

  GridContainers::Cartesian_gridpoint CreateCartesianGridpoint(
      const Eigen::Vector3d& atomA_pos,
      GridContainers::radial_grid& radial_grid,
      GridContainers::spherical_grid& spherical_grid, Index i_rad,
      Index i_sph) const;

  Eigen::VectorXd SSWpartition(const Eigen::VectorXd& rq_i,
                               const Eigen::MatrixXd& Rij) const;
  void SSWpartitionAtom(
      const QMMolecule& atoms,
      std::vector<GridContainers::Cartesian_gridpoint>& atomgrid, Index i_atom,
      const Eigen::MatrixXd& Rij) const;
  Eigen::MatrixXd CalcDistanceAtomsGridpoints(
      const QMMolecule& atoms,
      std::vector<GridContainers::Cartesian_gridpoint>& atomgrid) const;

  Index _totalgridsize;
  std::vector<GridBox> _grid_boxes;
  bool _density_set = false;
};

}  // namespace xtp
}  // namespace votca
#endif  // XTP_VXC_GRID_H
