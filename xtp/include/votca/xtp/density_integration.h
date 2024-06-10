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
#ifndef VOTCA_XTP_DENSITY_INTEGRATION_H
#define VOTCA_XTP_DENSITY_INTEGRATION_H

// Local VOTCA includes
#include "aobasis.h"
#include "eigen.h"
#include "regular_grid.h"
#include "vxc_grid.h"

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
  explicit DensityIntegration(const Grid& grid) : grid_(grid) {};

  double IntegrateDensity(const Eigen::MatrixXd& density_matrix);
  double IntegratePotential(const Eigen::Vector3d& rvector) const;
  Eigen::Vector3d IntegrateField(const Eigen::Vector3d& rvector) const;
  Eigen::MatrixXd IntegratePotential(const AOBasis& externalbasis) const;

  Gyrationtensor IntegrateGyrationTensor(const Eigen::MatrixXd& density_matrix);

  const std::vector<std::vector<double> >& getDensities() const {
    return densities_;
  }

 private:
  void SetupDensityContainer();
  const Grid grid_;

  std::vector<std::vector<double> > densities_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_DENSITY_INTEGRATION_H
