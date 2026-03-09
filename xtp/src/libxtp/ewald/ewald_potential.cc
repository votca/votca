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

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/ewald/ewald_potential.h"

namespace votca {
namespace xtp {
template <class Grid>
Ewald_Potential<Grid>::~Ewald_Potential() {}



template <class Grid>
Mat_p_Energy Ewald_Potential<Grid>::IntegrateEwald(Index basissize) const {

  Mat_p_Energy vewald = Mat_p_Energy(basissize,basissize);

#pragma omp parallel for schedule(guided) reduction(+ : vewald)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    const GridBox& box = grid_[i];
    if (!box.Matrixsize()) {
      continue;
    }
    double Eewald_box = 0.0;

    Eigen::MatrixXd Vewald_here =
        Eigen::MatrixXd::Zero(box.Matrixsize(),box.Matrixsize());
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    const std::vector<double>& pot = box.getPotentialValues();
    
    // iterate over gridpoints
    for (Index p = 0; p < box.size(); p++) {
      AOShell::AOValues ao = box.CalcAOValues(points[p]);
      const double weight = weights[p] * pot[p];
      if (weight < 1.e-20) {
        continue;  // skip the rest, if integrand is very small
      }

      Eewald_box += weight;
      Vewald_here.noalias() += weight * ao.values * ao.values.transpose();
    }
    box.AddtoBigMatrix(vewald.matrix(), Vewald_here);
    vewald.energy() += Eewald_box;
  }

  return Mat_p_Energy(vewald.energy(), vewald.matrix());
}

template class Ewald_Potential<Vxc_Grid>;

}  // namespace xtp
}  // namespace votca
