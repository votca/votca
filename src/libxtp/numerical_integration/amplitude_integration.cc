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

#include <votca/xtp/amplitude_integration.h>
#include <votca/xtp/aopotential.h>

namespace votca {
namespace xtp {

template <class Grid>
std::vector<std::vector<double> >
    AmplitudeIntegration<Grid>::IntegrateAmplitude(
        const Eigen::VectorXd& amplitude) {

  auto result = SetupAmplitudeContainer();

#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < _grid.getBoxesSize(); ++i) {
    const GridBox& box = _grid[i];
    const Eigen::VectorXd amplitude_here = box.ReadFromBigVector(amplitude);
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    // iterate over gridpoints
    for (Index p = 0; p < box.size(); p++) {
      Eigen::VectorXd ao = box.CalcAOValues(points[p]);
      result[i][p] = weights[p] * amplitude_here.dot(ao);
    }
  }
  return result;
}

template <class Grid>
std::vector<std::vector<double> >
    AmplitudeIntegration<Grid>::SetupAmplitudeContainer() {
  std::vector<std::vector<double> > amplitudes =
      std::vector<std::vector<double> >(_grid.getBoxesSize());
  for (Index i = 0; i < _grid.getBoxesSize(); i++) {
    amplitudes[i] = std::vector<double>(_grid[i].size(), 0.0);
  }
  return amplitudes;
}

template class AmplitudeIntegration<Regular_Grid>;

}  // namespace xtp
}  // namespace votca
