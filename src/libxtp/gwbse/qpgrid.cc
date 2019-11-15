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
  _roots.clear();
  PrintOptions();
  SetupGrid(frequencies);
  if (tools::globals::verbose) {
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " Evaluated sigma_c on QP grid" << std::flush;
  }
  FindRoots();
  PrintRoots();
  return CalculateQP_Energies(frequencies);
}

void QPGrid::SetupGrid(const Eigen::VectorXd frequencies) {
  const Eigen::VectorXd offset = Eigen::VectorXd::LinSpaced(_steps, -_range, +_range);
  for (Index i_node = 0; i_node < _steps; ++i_node) {
    Eigen::VectorXd _freqs_cur = frequencies.array() + offset[i_node];
    _freq.row(i_node) = _freqs_cur;
    _sigm.row(i_node) = _sigma.CalcCorrelationDiag(_freqs_cur);
  }
}

void QPGrid::FindRoots() {
  const Eigen::VectorXd temp = _sigma_x.diagonal() - _vxc.diagonal();
  for (Index i_qp = 0; i_qp < _qptotal; ++i_qp) {
    const double c = temp[i_qp] + _dft_energies[_opt.qpmin + i_qp];
    const Eigen::ArrayXd freq_cur = _freq.col(i_qp).array();
    const Eigen::ArrayXd sigm_cur = _sigm.col(i_qp).array();
    const Eigen::ArrayXd targ_cur = sigm_cur + c - freq_cur;
    for (Index i_node = 0; i_node < _steps - 1; ++i_node) {
      if (targ_cur[i_node] * targ_cur[i_node + 1] < 0.0) { // Sign change
        double dsdw = (sigm_cur[i_node + 1] - sigm_cur[i_node]) /
                      (freq_cur[i_node + 1] - freq_cur[i_node]);
        double dfdw = (targ_cur[i_node + 1] - targ_cur[i_node]) /
                      (freq_cur[i_node + 1] - freq_cur[i_node]);
        QPGrid::qproot root;
        root.i_qp = i_qp;
        root.i_node = i_node;
        root.value = freq_cur[i_node] - targ_cur[i_node] / dfdw; // Fixed-point estimate
        root.score = 1.0 / (1.0 - dsdw); // Pole weight Z
        _roots.push_back(root);
      }
    }
  }
}

Eigen::VectorXd QPGrid::CalculateQP_Energies(const Eigen::VectorXd frequencies) const {
  Eigen::VectorXd root_values = Eigen::VectorXd::Zero(_qptotal);
  Eigen::VectorXd root_scores = Eigen::VectorXd::Zero(_qptotal);
  for (Index i_root = 0; i_root < _roots.size(); ++i_root) {
    QPGrid::qproot root = _roots[i_root];
    if (root.score > 1e-5 && root.score > root_scores[root.i_qp]) {
      root_values[root.i_qp] = root.value;
      root_scores[root.i_qp] = root.score;
    }
  }
  return (root_scores.array() > 0).select(root_values, frequencies);
}

void QPGrid::PrintOptions() const {
  if (tools::globals::verbose) {
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " QP grid" << std::flush;
    XTP_LOG_SAVE(logINFO, _log)
        << boost::format(
            "Steps = %1$4d Range = %2$4d Delta = %3$+1.6f") %
            _steps % _range % _delta
        << std::flush;
  }
}

void QPGrid::PrintRoots() const {
  if (tools::globals::verbose) {
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " QP roots" << std::flush;
    for (Index i_root = 0; i_root < _roots.size(); ++i_root) {
      QPGrid::qproot root = _roots[i_root];
      if (root.i_qp <= _opt.homo) {
        XTP_LOG_SAVE(logINFO, _log)
            << boost::format(
                "Index = %1$4d Level = %2$4d Value = %3$+1.6f Ha Score = %4$+1.6f") %
                i_root % root.i_qp % root.value % root.score
            << std::flush;
      }
    }
  }
}

}  // namespace xtp
}  // namespace votca
