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

#include <votca/xtp/aopotential.h>
#include <votca/xtp/density_integration.h>
#include <votca/xtp/vxc_grid.h>

namespace votca {
namespace xtp {

template <class Grid>
double DensityIntegration<Grid>::IntegratePotential(
    const Eigen::Vector3d& rvector) const {

  double result = 0.0;
  assert(_densities.empty() && "Density not calculated");
  for (Index i = 0; i < _grid.getBoxesSize(); i++) {
    const std::vector<Eigen::Vector3d>& points = _grid[i].getGridPoints();
    const std::vector<double>& weights = _grid[i].getGridWeights();
    const std::vector<double>& densities = _densities[i];
    for (Index j = 0; j < _grid[i].size(); j++) {
      double charge = -weights[j] * densities[j];
      double dist = (points[j] - rvector).norm();
      result += charge / dist;
    }
  }
  return result;
}

template <class Grid>
Eigen::Vector3d DensityIntegration<Grid>::IntegrateField(
    const Eigen::Vector3d& rvector) const {

  Eigen::Vector3d result = Eigen::Vector3d::Zero();
  assert(_densities.empty() && "Density not calculated");
  for (Index i = 0; i < _grid.getBoxesSize(); i++) {
    const std::vector<Eigen::Vector3d>& points = _grid[i].getGridPoints();
    const std::vector<double>& weights = _grid[i].getGridWeights();
    const std::vector<double>& densities = _densities[i];
    for (Index j = 0; j < _grid[i].size(); j++) {
      double charge = -weights[j] * densities[j];
      Eigen::Vector3d r = points[j] - rvector;
      result += charge * r / std::pow(r.norm(), 3);  // x,y,z
    }
  }
  return result;
}

template <class Grid>
void DensityIntegration<Grid>::SetupDensityContainer() {
  _densities = std::vector<std::vector<double> >(_grid.getBoxesSize());
  for (Index i = 0; i < _grid.getBoxesSize(); i++) {
    _densities[i] = std::vector<double>(_grid[i].size(), 0.0);
  }
}

template <class Grid>
double DensityIntegration<Grid>::IntegrateDensity(
    const Eigen::MatrixXd& density_matrix) {

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<double> N_thread = std::vector<double>(nthreads, 0.0);
  SetupDensityContainer();

#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < _grid.getBoxesSize(); ++i) {
    const GridBox& box = _grid[i];
    double N_box = 0.0;
    const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(density_matrix);
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    // iterate over gridpoints
    for (Index p = 0; p < box.size(); p++) {
      Eigen::VectorXd ao = box.CalcAOValues(points[p]);
      double rho = (ao.transpose() * DMAT_here * ao)(0, 0);
      _densities[i][p] = rho;
      N_box += rho * weights[p];
    }
    N_thread[OPENMP::getThreadId()] += N_box;
  }
  double N = std::accumulate(N_thread.begin(), N_thread.end(), 0.0);
  return N;
}

template <class Grid>
Gyrationtensor DensityIntegration<Grid>::IntegrateGyrationTensor(
    const Eigen::MatrixXd& density_matrix) {

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<double> N_thread = std::vector<double>(nthreads, 0.0);
  std::vector<Eigen::Vector3d> centroid_thread =
      std::vector<Eigen::Vector3d>(nthreads, Eigen::Vector3d::Zero());
  std::vector<Eigen::Matrix3d> gyration_thread =
      std::vector<Eigen::Matrix3d>(nthreads, Eigen::Matrix3d::Zero());

  SetupDensityContainer();
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < _grid.getBoxesSize(); ++i) {
    const GridBox& box = _grid[i];
    double N_box = 0.0;
    Eigen::Vector3d centroid_box = Eigen::Vector3d::Zero();
    Eigen::Matrix3d gyration_box = Eigen::Matrix3d::Zero();
    const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(density_matrix);
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    // iterate over gridpoints
    for (Index p = 0; p < box.size(); p++) {
      Eigen::VectorXd ao = box.CalcAOValues(points[p]);
      double rho = (ao.transpose() * DMAT_here * ao).value();
      _densities[i][p] = rho;
      N_box += rho * weights[p];
      centroid_box += rho * weights[p] * points[p];
      gyration_box += rho * weights[p] * points[p] * points[p].transpose();
    }
    N_thread[OPENMP::getThreadId()] += N_box;
    centroid_thread[OPENMP::getThreadId()] += centroid_box;
    gyration_thread[OPENMP::getThreadId()] += gyration_box;
  }
  double N = std::accumulate(N_thread.begin(), N_thread.end(), 0.0);
  Eigen::Vector3d centroid =
      std::accumulate(centroid_thread.begin(), centroid_thread.end(),
                      Eigen::Vector3d::Zero().eval());
  Eigen::Matrix3d gyration =
      std::accumulate(gyration_thread.begin(), gyration_thread.end(),
                      Eigen::Matrix3d::Zero().eval());

  // Normalize
  centroid = centroid / N;
  gyration = gyration / N;
  gyration = gyration - centroid * centroid.transpose();
  Gyrationtensor gyro;
  gyro.mass = N;
  gyro.centroid = centroid;
  gyro.gyration = gyration;

  return gyro;
}
template <class Grid>
Eigen::MatrixXd DensityIntegration<Grid>::IntegratePotential(
    const AOBasis& externalbasis) const {
  Eigen::MatrixXd Potential = Eigen::MatrixXd::Zero(
      externalbasis.AOBasisSize(), externalbasis.AOBasisSize());

  assert(_densities.empty() && "Density not calculated");
  for (Index i = 0; i < _grid.getBoxesSize(); i++) {
    const std::vector<Eigen::Vector3d>& points = _grid[i].getGridPoints();
    const std::vector<double>& weights = _grid[i].getGridWeights();
    const std::vector<double>& densities = _densities[i];
    for (Index j = 0; j < _grid[i].size(); j++) {
      double weighteddensity = weights[j] * densities[j];
      if (weighteddensity < 1e-12) {
        continue;
      }
      AOMultipole esp;
      esp.FillPotential(externalbasis, points[j]);
      Potential += weighteddensity * esp.Matrix();
    }
  }
  return Potential;
}

template class DensityIntegration<Vxc_Grid>;

}  // namespace xtp
}  // namespace votca
