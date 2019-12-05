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
#include <votca/xtp/sphere_lebedev_rule.h>
#include <votca/xtp/vxc_grid.h>

namespace votca {
namespace xtp {

void Vxc_Grid::SortGridpointsintoBlocks(
    const std::vector<std::vector<GridContainers::Cartesian_gridpoint> >&
        grid) {
  const double boxsize = 1;  // 1 bohr

  Eigen::Array3d min =
      Eigen::Array3d::Ones() * std::numeric_limits<double>::max();
  Eigen::Array3d max =
      Eigen::Array3d::Ones() * std::numeric_limits<double>::min();

  for (const auto& atom_grid : grid) {
    for (const auto& gridpoint : atom_grid) {
      const Eigen::Vector3d& pos = gridpoint.grid_pos;
      max = max.max(pos.array()).eval();
      min = min.min(pos.array()).eval();
    }
  }

  Eigen::Array3d molextension = max - min;
  Eigen::Array<Index, 3, 1> numberofboxes =
      (molextension / boxsize).ceil().cast<Index>();

  std::vector<std::vector<
      std::vector<std::vector<const GridContainers::Cartesian_gridpoint*> > > >
      boxes;
  // creating temparray
  for (Index i = 0; i < numberofboxes.x(); i++) {
    std::vector<
        std::vector<std::vector<const GridContainers::Cartesian_gridpoint*> > >
        boxes_yz;
    for (Index j = 0; j < numberofboxes.y(); j++) {
      std::vector<std::vector<const GridContainers::Cartesian_gridpoint*> >
          boxes_z;
      for (Index k = 0; k < numberofboxes.z(); k++) {
        std::vector<const GridContainers::Cartesian_gridpoint*> box;
        box.reserve(100);
        boxes_z.push_back(box);
      }
      boxes_yz.push_back(boxes_z);
    }
    boxes.push_back(boxes_yz);
  }

  for (const auto& atomgrid : grid) {
    for (const auto& gridpoint : atomgrid) {
      Eigen::Vector3d pos = gridpoint.grid_pos - min.matrix();
      Eigen::Vector3d index = pos / boxsize;
      Index i_x = Index(index[0]);
      Index i_y = Index(index[1]);
      Index i_z = Index(index[2]);
      boxes[i_x][i_y][i_z].push_back(&gridpoint);
    }
  }
  for (auto& boxes_xy : boxes) {
    for (auto& boxes_z : boxes_xy) {
      for (auto& box : boxes_z) {
        if (box.empty()) {
          continue;
        }
        GridBox gridbox;
        for (const auto& point : box) {
          gridbox.addGridPoint(*point);
        }
        _grid_boxes.push_back(gridbox);
      }
    }
  }
  return;
}

void Vxc_Grid::FindSignificantShells(const AOBasis& basis) {

#pragma omp parallel for
  for (Index i = 0; i < getBoxesSize(); i++) {
    _grid_boxes[i].FindSignificantShells(basis);
  }

  std::vector<GridBox> grid_boxes_copy;
  // use vector of bool to indicate if a gridbox has already been merged into
  // another
  std::vector<bool> Merged = std::vector<bool>(_grid_boxes.size(), false);
  for (Index i = 0; i < Index(_grid_boxes.size()); i++) {
    if (Merged[i]) {
      continue;
    }
    GridBox box = _grid_boxes[i];
    if (box.Shellsize() < 1) {
      continue;
    }
    Merged[i] = true;
    for (Index j = i + 1; j < Index(_grid_boxes.size()); j++) {
      if (GridBox::compareGridboxes(_grid_boxes[i], _grid_boxes[j])) {
        Merged[j] = true;
        box.addGridBox(_grid_boxes[j]);
      }
    }
    grid_boxes_copy.push_back(box);
  }

  _totalgridsize = 0;
  for (auto& box : grid_boxes_copy) {
    _totalgridsize += box.size();
    box.PrepareForIntegration();
  }
  _grid_boxes = grid_boxes_copy;
}

std::vector<const Eigen::Vector3d*> Vxc_Grid::getGridpoints() const {
  std::vector<const Eigen::Vector3d*> gridpoints;
  gridpoints.reserve(this->getGridSize());
  for (const auto& box : _grid_boxes) {
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    for (const Eigen::Vector3d& point : points) {
      gridpoints.push_back(&point);
    }
  }
  return gridpoints;
}

Eigen::MatrixXd Vxc_Grid::CalcInverseAtomDist(const QMMolecule& atoms) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(atoms.size(), atoms.size());
#pragma omp parallel for
  for (Index i = 0; i < atoms.size(); ++i) {
    const Eigen::Vector3d& pos_a = atoms[i].getPos();
    for (Index j = 0; j < i; ++j) {
      const Eigen::Vector3d& pos_b = atoms[j].getPos();
      result(j, i) = 1 / (pos_a - pos_b).norm();
    }
  }
  return result + result.transpose();
}

Index Vxc_Grid::UpdateOrder(LebedevGrid& sphericalgridofElement, Index maxorder,
                            std::vector<double>& PruningIntervals,
                            double r) const {
  Index order;
  Index maxindex = sphericalgridofElement.getIndexFromOrder(maxorder);
  if (maxindex == 1) {
    // smallest possible grid anyway, nothing to do
    order = maxorder;
  } else if (maxindex == 2) {
    // only three intervals
    if (r < PruningIntervals[0]) {
      order = sphericalgridofElement.getOrderFromIndex(1);  // 1;
    } else if ((r >= PruningIntervals[0]) && (r < PruningIntervals[3])) {
      order = sphericalgridofElement.getOrderFromIndex(2);
    } else {
      order = sphericalgridofElement.getOrderFromIndex(1);
    }  // maxorder == 2
  } else {
    // five intervals
    if (r < PruningIntervals[0]) {
      order = sphericalgridofElement.getOrderFromIndex(2);
    } else if ((r >= PruningIntervals[0]) && (r < PruningIntervals[1])) {
      order = sphericalgridofElement.getOrderFromIndex(4);
    } else if ((r >= PruningIntervals[1]) && (r < PruningIntervals[2])) {
      order =
          sphericalgridofElement.getOrderFromIndex(std::max(maxindex - 1, 4l));
    } else if ((r >= PruningIntervals[2]) && (r < PruningIntervals[3])) {
      order = maxorder;
    } else {
      order =
          sphericalgridofElement.getOrderFromIndex(std::max(maxindex - 1, 1l));
    }
  }
  return order;
}

GridContainers::Cartesian_gridpoint Vxc_Grid::CreateCartesianGridpoint(
    const Eigen::Vector3d& atomA_pos, GridContainers::radial_grid& radial_grid,
    GridContainers::spherical_grid& spherical_grid, Index i_rad,
    Index i_sph) const {
  GridContainers::Cartesian_gridpoint gridpoint;
  double p = spherical_grid.phi[i_sph];
  double t = spherical_grid.theta[i_sph];
  const Eigen::Vector3d s =
      Eigen::Vector3d{sin(p) * cos(t), sin(p) * sin(t), cos(p)};
  double r = radial_grid.radius[i_rad];
  gridpoint.grid_pos = atomA_pos + r * s;
  gridpoint.grid_weight =
      radial_grid.weight[i_rad] * spherical_grid.weight[i_sph];
  return gridpoint;
}

Eigen::MatrixXd Vxc_Grid::CalcDistanceAtomsGridpoints(
    const QMMolecule& atoms,
    std::vector<GridContainers::Cartesian_gridpoint>& atomgrid) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(atoms.size(), atomgrid.size());
#pragma omp parallel for
  for (Index i = 0; i < atoms.size(); ++i) {
    const Eigen::Vector3d& atom_pos = atoms[i].getPos();
    for (Index j = 0; j < Index(atomgrid.size()); ++j) {
      const auto& gridpoint = atomgrid[j];
      result(i, j) = (atom_pos - gridpoint.grid_pos).norm();
    }
  }
  return result;
}

void Vxc_Grid::SSWpartitionAtom(
    const QMMolecule& atoms,
    std::vector<GridContainers::Cartesian_gridpoint>& atomgrid, Index i_atom,
    const Eigen::MatrixXd& Rij) const {
  Eigen::MatrixXd AtomGridDist = CalcDistanceAtomsGridpoints(atoms, atomgrid);

#pragma omp parallel for schedule(guided)
  for (Index i_grid = 0; i_grid < Index(atomgrid.size()); i_grid++) {
    Eigen::VectorXd p = SSWpartition(AtomGridDist.col(i_grid), Rij);
    // check weight sum
    double wsum = p.sum();
    if (wsum != 0.0) {
      // update the weight of this grid point
      atomgrid[i_grid].grid_weight *= p[i_atom] / wsum;
    } else {
      std::cerr << "\nSum of partition weights of grid point " << i_grid
                << " of atom " << i_atom << " is zero! ";
      throw std::runtime_error("\nThis should never happen!");
    }
  }  // partition weight for each gridpoint
}

void Vxc_Grid::GridSetup(const std::string& type, const QMMolecule& atoms,
                         const AOBasis& basis) {
  GridContainers initialgrids;
  // get radial grid per element
  EulerMaclaurinGrid radialgridofElement;
  initialgrids.radial_grids = radialgridofElement.CalculateAtomicRadialGrids(
      basis, atoms, type);  // this checks out 1:1 with NWChem results! AWESOME
  LebedevGrid sphericalgridofElement;
  initialgrids.spherical_grids =
      sphericalgridofElement.CalculateSphericalGrids(atoms, type);

  // for the partitioning, we need all inter-center distances later, stored in
  // matrix
  Eigen::MatrixXd Rij = CalcInverseAtomDist(atoms);
  std::vector<std::vector<GridContainers::Cartesian_gridpoint> > grid;

  for (Index i_atom = 0; i_atom < atoms.size(); ++i_atom) {
    const QMAtom& atom = atoms[i_atom];

    const Eigen::Vector3d& atomA_pos = atom.getPos();
    const std::string& name = atom.getElement();
    GridContainers::radial_grid radial_grid =
        initialgrids.radial_grids.at(name);
    GridContainers::spherical_grid spherical_grid =
        initialgrids.spherical_grids.at(name);

    // maximum order (= number of points) in spherical integration grid
    Index maxorder = sphericalgridofElement.Type2MaxOrder(name, type);
    // for pruning of integration grid, get interval boundaries for this element
    std::vector<double> PruningIntervals =
        radialgridofElement.CalculatePruningIntervals(name);
    Index current_order = 0;
    // for each radial value
    std::vector<GridContainers::Cartesian_gridpoint> atomgrid;
    for (Index i_rad = 0; i_rad < radial_grid.radius.size(); i_rad++) {
      double r = radial_grid.radius[i_rad];

      // which Lebedev order for this point?
      Index order =
          UpdateOrder(sphericalgridofElement, maxorder, PruningIntervals, r);
      // get new spherical grid, if order changed
      if (order != current_order) {
        spherical_grid = sphericalgridofElement.CalculateUnitSphereGrid(order);
        current_order = order;
      }

      for (Index i_sph = 0; i_sph < spherical_grid.phi.size(); i_sph++) {
        GridContainers::Cartesian_gridpoint gridpoint =
            CreateCartesianGridpoint(atomA_pos, radial_grid, spherical_grid,
                                     i_rad, i_sph);
        atomgrid.push_back(gridpoint);
      }  // spherical gridpoints
    }    // radial gridpoint

    SSWpartitionAtom(atoms, atomgrid, i_atom, Rij);
    // now remove points from the grid with negligible weights
    std::vector<GridContainers::Cartesian_gridpoint> atomgrid_cleanedup;
    for (const auto& point : atomgrid) {
      if (point.grid_weight > 1e-13) {
        atomgrid_cleanedup.push_back(point);
      }
    }
    grid.push_back(atomgrid_cleanedup);
  }  // atoms
  SortGridpointsintoBlocks(grid);
  FindSignificantShells(basis);
  return;
}

Eigen::VectorXd Vxc_Grid::SSWpartition(const Eigen::VectorXd& rq_i,
                                       const Eigen::MatrixXd& Rij) const {
  const double ass = 0.725;
  // initialize partition vector to 1.0
  Eigen::VectorXd p = Eigen::VectorXd::Ones(rq_i.size());
  const double tol_scr = 1e-10;
  const double leps = 1e-6;
  // go through centers
  for (Index i = 1; i < rq_i.size(); i++) {
    double rag = rq_i(i);
    // through all other centers (one-directional)
    for (Index j = 0; j < i; j++) {
      if ((std::abs(p[i]) > tol_scr) || (std::abs(p[j]) > tol_scr)) {
        double mu = (rag - rq_i(j)) * Rij(j, i);
        if (mu > ass) {
          p[i] = 0.0;
        } else if (mu < -ass) {
          p[j] = 0.0;
        } else {
          double sk;
          if (std::abs(mu) < leps) {
            sk = -1.88603178008 * mu + 0.5;
          } else {
            sk = erf1c(mu);
          }
          if (mu > 0.0) {
            sk = 1.0 - sk;
          }
          p[j] = p[j] * sk;
          p[i] = p[i] * (1.0 - sk);
        }
      }
    }
  }
  return p;
}

double Vxc_Grid::erf1c(double x) const {
  const static double alpha_erf1 = 1.0 / 0.30;
  return 0.5 * std::erfc(std::abs(x / (1.0 - x * x)) * alpha_erf1);
}

}  // namespace xtp
}  // namespace votca
