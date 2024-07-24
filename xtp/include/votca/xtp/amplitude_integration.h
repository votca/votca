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
#ifndef VOTCA_XTP_AMPLITUDE_INTEGRATION_H
#define VOTCA_XTP_AMPLITUDE_INTEGRATION_H

// Local VOTCA includes
#include "aobasis.h"
#include "eigen.h"
#include "regular_grid.h"

namespace votca {
namespace xtp {

template <class Grid>
class AmplitudeIntegration {
 public:
  explicit AmplitudeIntegration(const Grid& grid) : grid_(grid) {};

  std::vector<std::vector<double> > IntegrateAmplitude(
      const Eigen::VectorXd& amplitude);

 private:
  std::vector<std::vector<double> > SetupAmplitudeContainer();
  const Grid grid_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_AMPLITUDE_INTEGRATION_H
