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
#ifndef XTP_DENSITY_INTEGRATION_H
#define XTP_DENSITY_INTEGRATION_H

#include <votca/xtp/aobasis.h>
#include <votca/xtp/eigen.h>
namespace votca {
namespace xtp {

struct Gyrationtensor {
  double mass;
  Eigen::Vector3d centroid;
  Eigen::Matrix3d gyration;
};

template <class Grid>
class DensityIntegration {
 public:
  DensityIntegration(Grid grid) : _grid(grid){};

  double IntegrateDensity(const Eigen::MatrixXd& density_matrix);
  double IntegratePotential(const Eigen::Vector3d& rvector) const;
  Eigen::Vector3d IntegrateField(const Eigen::Vector3d& rvector) const;
  Eigen::MatrixXd IntegratePotential(const AOBasis& externalbasis) const;

  Gyrationtensor IntegrateGyrationTensor(const Eigen::MatrixXd& density_matrix);

 private:
  void SetupDensityContainer();
  const Grid _grid;

  std::vector<std::vector<double> > _densities;
};

}  // namespace xtp
}  // namespace votca
#endif  // XTP_NUMERICAL_INTEGRATION_H
