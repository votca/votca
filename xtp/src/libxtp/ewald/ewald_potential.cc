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
#include "votca/xtp/ewald_potential.h"

namespace votca {
namespace xtp {
template <class Grid>
Ewald_Potential<Grid>::~Ewald_Potential() {}



template <class Grid>
Mat_p_Energy Ewald_Potential<Grid>::IntegrateVXC(
    const Eigen::MatrixXd& density_matrix) const {

  assert(density_matrix.isApprox(density_matrix.transpose()) &&
         "Density matrix has to be symmetric!");
  Mat_p_Energy vxc = Mat_p_Energy(density_matrix.rows(), density_matrix.cols());

#pragma omp parallel for schedule(guided) reduction(+ : vxc)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    const GridBox& box = grid_[i];
    if (!box.Matrixsize()) {
      continue;
    }
    double EXC_box = 0.0;
    // two because we have to use the density matrix and its transpose
    const Eigen::MatrixXd DMAT_here = 2 * box.ReadFromBigMatrix(density_matrix);
    double cutoff =
        1.e-40 / double(density_matrix.rows()) / double(density_matrix.rows());
    if (DMAT_here.cwiseAbs2().maxCoeff() < cutoff) {
      continue;
    }
    Eigen::MatrixXd Vxc_here =
        Eigen::MatrixXd::Zero(DMAT_here.rows(), DMAT_here.cols());
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();

    // iterate over gridpoints
    for (Index p = 0; p < box.size(); p++) {
      AOShell::AOValues ao = box.CalcAOValues(points[p]);
      Eigen::VectorXd temp = ao.values.transpose() * DMAT_here;
      double rho = 0.5 * temp.dot(ao.values);
      const double weight = weights[p];
      if (rho * weight < 1.e-20) {
        continue;  // skip the rest, if density is very small
      }
      const Eigen::Vector3d rho_grad = temp.transpose() * ao.derivatives;

      typename Vxc_Potential<Grid>::XC_entry xc =
          EvaluateXC(rho, rho_grad.squaredNorm());
      EXC_box += weight * rho * xc.f_xc;
      auto grad = ao.derivatives * rho_grad;
      temp.noalias() =
          weight * (0.5 * xc.df_drho * ao.values + 2.0 * xc.df_dsigma * grad);
      Vxc_here.noalias() += temp * ao.values.transpose();
    }
    box.AddtoBigMatrix(vxc.matrix(), Vxc_here);
    vxc.energy() += EXC_box;
  }

  return Mat_p_Energy(vxc.energy(), vxc.matrix() + vxc.matrix().transpose());
}

template class Ewald_Potential<Vxc_Grid>;

}  // namespace xtp
}  // namespace votca
