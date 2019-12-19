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
#ifndef XTP_AMPLITUDE_INTEGRATION_H
#define XTP_AMPLITUDE_INTEGRATION_H

#include <votca/xtp/aobasis.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/regular_grid.h>
namespace votca {
namespace xtp {

template <class Grid>
class AmplitudeIntegration {
 public:
  explicit AmplitudeIntegration(const Grid& grid) : _grid(grid){};

  std::vector<std::vector<double> > IntegrateAmplitude(
      const Eigen::VectorXd& amplitude);

 private:
  std::vector<std::vector<double> > SetupAmplitudeContainer();
  const Grid _grid;
};

}  // namespace xtp
}  // namespace votca
#endif  // XTP_NUMERICAL_INTEGRATION_H
