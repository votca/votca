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

namespace votca {
namespace xtp {

Eigen::VectorXd QPGrid::Evaluate(const Eigen::VectorXd frequencies) {
  const Index qptotal = _opt.qpmax - _opt.qpmin + 1;
  const double range = _opt.spacing * (double)(_opt.steps - 1) / 2.0;
  Eigen::VectorXd offset = Eigen::VectorXd::LinSpaced(_opt.steps, -range, range);
  Eigen::MatrixXd sigc_mat = Eigen::MatrixXd::Zero(qptotal, _opt.steps);
  for (Index i_node = 0; i_node < _opt.steps; ++i_node) {
    sigc_mat.col(i_node) = _sigma.CalcCorrelationDiag(frequencies.array() + offset[i_node]);
  }
  Eigen::VectorXd roots = frequencies;
  for (Index level = 0; level < qptotal; ++level) {
    Eigen::VectorXd freq_cur = frequencies[level] + offset.array();
    Eigen::VectorXd sigc_cur = sigc_mat.row(level);
    Eigen::VectorXd targ_cur = sigc_cur.array() + _qpOffset[level] - freq_cur.array();
    double root_opt = 0.0;
    double pole_weight_max = -1.0;
    for (Index i_node = 0; i_node < _opt.steps - 1; ++i_node) {
      if (targ_cur[i_node] * targ_cur[i_node + 1] < 0.0) {  // Sign change
        double dsigc_dfreq = (sigc_cur[i_node + 1] - sigc_cur[i_node]) /
                             (freq_cur[i_node + 1] - freq_cur[i_node]);
        double dtarg_dfreq = (targ_cur[i_node + 1] - targ_cur[i_node]) /
                             (freq_cur[i_node + 1] - freq_cur[i_node]);
        double root = freq_cur[i_node] - targ_cur[i_node] / dtarg_dfreq;  // Fixed-point estimate
        double pole_weight = 1.0 / (1.0 - dsigc_dfreq);  // Pole weight Z in (0, 1)
        if (pole_weight >= 1e-5 && pole_weight > pole_weight_max) {
          root_opt = root;
          pole_weight_max = pole_weight;
        }
      }
    }
    if (pole_weight_max >= 0.0) {
      roots[level] = root_opt;
    }
  }
  return roots;
}

}  // namespace xtp
}  // namespace votca
