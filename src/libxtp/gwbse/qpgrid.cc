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

#include "votca/xtp/qpgrid.h"
#include <boost/format.hpp>
#include <votca/tools/globals.h>

namespace votca {
namespace xtp {

Eigen::VectorXd QPGrid::Evaluate(const Eigen::VectorXd frequencies) {
  const Eigen::VectorXd sigx_min_vxc = _sigma_x.diagonal() - _vxc.diagonal();
  const Index qptotal = _opt.qpmax - _opt.qpmin + 1;
  const double delta = (2.0 * _opt.range) / (_opt.steps - 1.0);
  Eigen::MatrixXd sigmc_mat = Eigen::MatrixXd::Zero(qptotal, _opt.steps);
  for (Index i_node = 0; i_node < _opt.steps; ++i_node) {
    sigmc_mat.col(i_node) = _sigma.CalcCorrelationDiag(frequencies.array() - _opt.range + i_node * delta);
  }
  Eigen::VectorXd roots = frequencies;
  for (Index level = 0; level < qptotal; ++level) {
    const double freq_cur = frequencies[level];
    Eigen::VectorXd grid_freq = freq_cur + Eigen::VectorXd::LinSpaced(_opt.steps, -_opt.range, _opt.range).array();
    Eigen::VectorXd grid_sigc = sigmc_mat.row(level);
    Eigen::VectorXd grid_targ = grid_sigc.array() + sigx_min_vxc[level] + _dft_energies[_opt.qpmin + level] - freq_cur;
    double root_best = 0.0;
    double weight_best = -1.0;
    for (Index i_node = 0; i_node < _opt.steps - 1; ++i_node) {
      if (grid_targ[i_node] * grid_targ[i_node + 1] < 0.0) { // Sign change
        double dsigc_dw = (grid_sigc[i_node + 1] - grid_sigc[i_node]) /
                          (grid_freq[i_node + 1] - grid_freq[i_node]);
        double dtarg_dw = (grid_targ[i_node + 1] - grid_targ[i_node]) /
                          (grid_freq[i_node + 1] - grid_freq[i_node]);
        double root = grid_freq[i_node] - grid_targ[i_node] / dtarg_dw; // Fixed-point estimate
        double weight = 1.0 / (1.0 - dsigc_dw); // Pole weight Z in (0, 1)
        if (weight >= 1e-5 && weight > weight_best) {
          root_best = root;
          weight_best = weight;
        }
      }
    }
    if (weight_best >= 0.0) {
      roots[level] = root_best;
    }
  }
  return roots;
}

}  // namespace xtp
}  // namespace votca
