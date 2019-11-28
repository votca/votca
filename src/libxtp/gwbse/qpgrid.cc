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

Eigen::VectorXd QPGrid::Evaluate(const Eigen::VectorXd& gwa_energies) {
  const Index qptotal = _opt.qpmax - _opt.qpmin + 1;
  const double range = _opt.spacing * (double)(_opt.steps - 1) / 2.0;
  Eigen::VectorXd offset = Eigen::VectorXd::LinSpaced(_opt.steps, -range, range);
  Eigen::MatrixXd sigc_mat = Eigen::MatrixXd::Zero(qptotal, _opt.steps);
  for (Index i_node = 0; i_node < _opt.steps; ++i_node) {
    sigc_mat.col(i_node) = _sigma.CalcCorrelationDiag(gwa_energies.array() + offset[i_node]);
  }
  Eigen::VectorXd gwa_energies_new = gwa_energies;
  for (Index level = 0; level < qptotal; ++level) {
    const Eigen::VectorXd freq = gwa_energies[level] + offset.array();
    const Eigen::VectorXd sigc = sigc_mat.row(level);
    const Eigen::VectorXd targ = sigc.array() + _energy_intercept[level] - freq.array();
    double qp_energy;
    bool found = FindQPEnergy(freq, sigc, targ, qp_energy);
    if (found) {
      gwa_energies_new[level] = qp_energy;
    }
  }
  return gwa_energies_new;
}

bool QPGrid::FindQPEnergy(const Eigen::VectorXd& freq, const Eigen::VectorXd& sigc, const Eigen::VectorXd& targ, double& qp_energy) const {
  double pole_weight_max = -1.0;
  for (Index i_node = 0; i_node < _opt.steps - 1; ++i_node) {
    if (targ[i_node] * targ[i_node + 1] < 0.0) {  // Sign change
      double dsigc_dfreq = (sigc[i_node + 1] - sigc[i_node]) /
                           (freq[i_node + 1] - freq[i_node]);
      double dtarg_dfreq = (targ[i_node + 1] - targ[i_node]) /
                           (freq[i_node + 1] - freq[i_node]);
      double pole = freq[i_node] - targ[i_node] / dtarg_dfreq;  // Fixed-point estimate
      double pole_weight = 1.0 / (1.0 - dsigc_dfreq);  // Pole weight Z in (0, 1)
      if (pole_weight >= 1e-5 && pole_weight > pole_weight_max) {
        qp_energy = pole;
        pole_weight_max = pole_weight;
      }
    }
  }
  return pole_weight_max >= 0.0;
}

}  // namespace xtp
}  // namespace votca
