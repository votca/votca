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
#include <iostream>
#include <numeric>
#include <sstream>
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

Eigen::MatrixXd RPA_RI_Reduced::NormalizeMetric(const Eigen::MatrixXd& M) {
  const double fnorm = M.norm();
  if (fnorm < 1e-16) {
    return M;
  }
  return M / fnorm;
}

Eigen::MatrixXd RPA_RI_Reduced::BuildStaticOrDynamicCompressionMatrix() const {
  const std::vector<double> grid = ImagFrequencyGrid();

  if (!opt_red_.dynamic_basis || grid.empty()) {
    return BuildPiImag(0.0);
  }

  const Index n_aux = Mmn_.auxsize();
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n_aux, n_aux);

  // Symmetric positive semidefinite dynamical metric:
  //   C = sum_k Pi(iw_k)^2
  for (double omega : grid) {
    const Eigen::MatrixXd Piw = BuildPiImag(omega);
    C.noalias() += Piw * Piw;
  }

  if (C.norm() < 1e-16) {
    return BuildPiImag(0.0);
  }

  return C;
}

Eigen::MatrixXd RPA_RI_Reduced::BuildSigmaAwareCompressionMatrix() const {
  const RPAWindow w = GetWindow();
  const Index n_aux = Mmn_.auxsize();
  const Index qptotal = (qpmax_ >= qpmin_) ? (qpmax_ - qpmin_ + 1) : 0;
  const Index rpatotal = w.n_occ + w.n_unocc;

  Eigen::MatrixXd Csigma = Eigen::MatrixXd::Zero(n_aux, n_aux);

  if (qptotal <= 0 || rpatotal <= 0) {
    return Csigma;
  }

#pragma omp parallel
  {
    Eigen::MatrixXd Cthread = Eigen::MatrixXd::Zero(n_aux, n_aux);

#pragma omp for schedule(dynamic)
    for (Index iq = 0; iq < qptotal; ++iq) {
      const Index level = qpmin_ + iq;
      const Eigen::MatrixXd& Mim = Mmn_[level - rpamin_];

      for (Index m = 0; m < rpatotal; ++m) {
        const Eigen::VectorXd c = Mim.row(m).transpose();
        Cthread.noalias() += c * c.transpose();
      }
    }

#pragma omp critical
    Csigma += Cthread;
  }

  return Csigma;
}

Eigen::MatrixXd RPA_RI_Reduced::BuildCompressionMatrix() const {
  Eigen::MatrixXd Cbase = BuildStaticOrDynamicCompressionMatrix();

  if (!opt_red_.sigma_aware_basis) {
    return Cbase;
  }

  Eigen::MatrixXd Csigma = BuildSigmaAwareCompressionMatrix();
  if (Csigma.norm() < 1e-16) {
    return Cbase;
  }

  const double mix = std::max(0.0, std::min(1.0, opt_red_.sigma_mix));

  if (opt_red_.normalize_metric_components) {
    Cbase = NormalizeMetric(Cbase);
    Csigma = NormalizeMetric(Csigma);
  }

  Eigen::MatrixXd C = (1.0 - mix) * Cbase + mix * Csigma;

  if (C.norm() < 1e-16) {
    return BuildStaticOrDynamicCompressionMatrix();
  }

  return C;
}

void RPA_RI_Reduced::BuildReducedBasis() {
  const Eigen::MatrixXd compression = BuildCompressionMatrix();

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(compression);
  if (es.info() != Eigen::Success) {
    throw std::runtime_error(
        "RPA_RI_Reduced::BuildReducedBasis failed to diagonalize compression matrix.");
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
    keep.push_back(order.front());
  }

  U_.resize(compression.rows(), static_cast<Index>(keep.size()));
  for (Index i = 0; i < static_cast<Index>(keep.size()); ++i) {
    U_.col(i) = evecs.col(keep[static_cast<std::size_t>(i)]);
  }

  if (opt_red_.print_basis_diagnostics) {
    std::ostringstream oss;
    oss << "[RI-REDUCED BASIS] mode = "
        << (opt_red_.dynamic_basis ? "dynamic" : "static")
        << "  sigma_aware = " << (opt_red_.sigma_aware_basis ? "true" : "false")
        << "  sigma_mix = " << opt_red_.sigma_mix
        << "  rank = " << U_.cols()
        << "  imag_omega_points = " << opt_red_.imag_omega_points
        << "  imag_omega_max = " << opt_red_.imag_omega_max << "\n";

    const std::vector<double> grid = ImagFrequencyGrid();
    for (double omega : grid) {
      const auto diag = CompareWcImag(omega);
      oss << "  omega = " << diag.omega
          << "  ||Wc_full||_F = " << diag.norm_full
          << "  ||Wc_proj||_F = " << diag.norm_proj
          << "  ||diff||_F = " << diag.abs_diff
          << "  rel = " << diag.rel_diff << "\n";
    }

    std::cout << oss.str() << std::flush;
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

  return (I - pi_red).inverse() - I;
}

Eigen::MatrixXd RPA_RI_Reduced::BuildWcImag(double omega) const {
  const Eigen::MatrixXd pi = BuildPiImag(omega);
  const Index dim = pi.rows();
  const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim, dim);

  return (I - pi).inverse() - I;
}

Eigen::MatrixXd RPA_RI_Reduced::BuildProjectedReducedWcImag(double omega) const {
  if (U_.cols() == 0) {
    throw std::runtime_error(
        "RPA_RI_Reduced::BuildProjectedReducedWcImag called before BuildReducedBasis.");
  }

  const Eigen::MatrixXd wc_red = BuildReducedWcImag(omega);
  return U_ * wc_red * U_.transpose();
}

RPA_RI_Reduced::wc_diagnostic RPA_RI_Reduced::CompareWcImag(double omega) const {
  const Eigen::MatrixXd wc_full = BuildWcImag(omega);
  const Eigen::MatrixXd wc_proj = BuildProjectedReducedWcImag(omega);
  const Eigen::MatrixXd diff = wc_full - wc_proj;

  wc_diagnostic diag;
  diag.omega = omega;
  diag.norm_full = wc_full.norm();
  diag.norm_proj = wc_proj.norm();
  diag.abs_diff = diff.norm();
  diag.rel_diff =
      (diag.norm_full > 1e-14) ? (diag.abs_diff / diag.norm_full) : diag.abs_diff;

  return diag;
}

}  // namespace xtp
}  // namespace votca