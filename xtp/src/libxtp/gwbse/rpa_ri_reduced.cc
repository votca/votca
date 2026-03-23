/*
 *            Copyright 2009-2026 The VOTCA Development Team
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

#include "votca/xtp/rpa_ri_reduced.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

RPA_RI_Reduced::RPAWindow RPA_RI_Reduced::GetWindow() const {
  const Index lumo = homo_ + 1;
  const Index n_occ = lumo - rpamin_;
  const Index n_unocc = rpamax_ - homo_;
  return {lumo, n_occ, n_unocc};
}

std::vector<double> RPA_RI_Reduced::ImagFrequencyGrid() const {
  std::vector<double> grid;

  if (opt_red_.imag_omega_points <= 0) {
    return grid;
  }

  if (opt_red_.imag_omega_points == 1) {
    grid.push_back(0.0);
    return grid;
  }

  grid.reserve(static_cast<std::size_t>(opt_red_.imag_omega_points));
  const double wmax = opt_red_.imag_omega_max;
  const double denom = static_cast<double>(opt_red_.imag_omega_points - 1);

  for (Index i = 0; i < opt_red_.imag_omega_points; ++i) {
    grid.push_back(wmax * static_cast<double>(i) / denom);
  }

  return grid;
}

Eigen::MatrixXd RPA_RI_Reduced::BuildPiImag(double omega) const {
  if (energies_.size() == 0) {
    throw std::runtime_error(
        "RPA_RI_Reduced::BuildPiImag called without RPA input energies.");
  }

  const RPAWindow w = GetWindow();
  const Index n_aux = Mmn_.auxsize();

  Eigen::MatrixXd Pi = Eigen::MatrixXd::Zero(n_aux, n_aux);

#pragma omp parallel
  {
    Eigen::MatrixXd Pi_thread = Eigen::MatrixXd::Zero(n_aux, n_aux);

#pragma omp for schedule(dynamic)
    for (Index v = 0; v < w.n_occ; ++v) {
      const double eps_v = energies_(v);

      // rows correspond to n in [rpamin, rpamax]
      // middleRows(n_occ, n_unocc) == virtual block c for fixed occupied v
      const Eigen::MatrixXd Mvc = Mmn_[v].middleRows(w.n_occ, w.n_unocc);

      Eigen::VectorXd weight(w.n_unocc);
      for (Index c = 0; c < w.n_unocc; ++c) {
        const double delta = energies_(w.n_occ + c) - eps_v;
        weight(c) = -4.0 * delta / (omega * omega + delta * delta);
      }

      Pi_thread.noalias() += Mvc.transpose() * weight.asDiagonal() * Mvc;
    }

#pragma omp critical
    Pi += Pi_thread;
  }

  return Pi;
}

void RPA_RI_Reduced::BuildReducedBasis() {
  const Eigen::MatrixXd Pi0 = BuildPiImag(0.0);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Pi0);
  if (es.info() != Eigen::Success) {
    throw std::runtime_error(
        "RPA_RI_Reduced::BuildReducedBasis failed to diagonalize Pi(0).");
  }

  const Eigen::VectorXd evals = es.eigenvalues();
  const Eigen::MatrixXd evecs = es.eigenvectors();

  std::vector<Index> order(static_cast<std::size_t>(evals.size()));
  std::iota(order.begin(), order.end(), 0);

  std::sort(order.begin(), order.end(),
            [&](Index a, Index b) { return std::abs(evals(a)) > std::abs(evals(b)); });

  std::vector<Index> keep;
  keep.reserve(order.size());

  for (Index idx : order) {
    if (std::abs(evals(idx)) < opt_red_.basis_threshold) {
      continue;
    }
    keep.push_back(idx);
    if (opt_red_.max_rank > 0 &&
        static_cast<Index>(keep.size()) >= opt_red_.max_rank) {
      break;
    }
  }

  if (keep.empty()) {
    // Keep at least the dominant mode so the reduced space is never empty.
    keep.push_back(order.front());
  }

  U_.resize(Pi0.rows(), static_cast<Index>(keep.size()));
  for (Index i = 0; i < static_cast<Index>(keep.size()); ++i) {
    U_.col(i) = evecs.col(keep[static_cast<std::size_t>(i)]);
  }
}

Eigen::MatrixXd RPA_RI_Reduced::BuildReducedPiImag(double omega) const {
  if (U_.cols() == 0) {
    throw std::runtime_error(
        "RPA_RI_Reduced::BuildReducedPiImag called before BuildReducedBasis.");
  }

  return U_.transpose() * BuildPiImag(omega) * U_;
}

Eigen::MatrixXd RPA_RI_Reduced::BuildReducedWcImag(double omega) const {
  const Eigen::MatrixXd pi_red = BuildReducedPiImag(omega);
  const Index dim = pi_red.rows();
  const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim, dim);

  // Correlation part of reduced screened interaction in the reduced basis
  return (I - pi_red).inverse() - I;
}

}  // namespace xtp
}  // namespace votca