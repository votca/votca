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

#pragma once
#ifndef VOTCA_XTP_RPA_RI_REDUCED_H
#define VOTCA_XTP_RPA_RI_REDUCED_H

#include <vector>

#include "eigen.h"

namespace votca {
namespace xtp {

class TCMatrix_gwbse;

/**
 * Reduced RI-space RPA screening model.
 *
 * This class avoids the explicit diagonalization of the full transition-space
 * RPA Hamiltonian. Instead, it constructs the independent-particle
 * polarizability directly in the Coulomb-whitened auxiliary RI basis already
 * stored in TCMatrix_gwbse, and compresses that screening response into a
 * truncated auxiliary subspace U.
 */
class RPA_RI_Reduced {
 public:
  struct options {
    double imag_omega_max = 2.0;
    Index imag_omega_points = 6;
    double basis_threshold = 1e-8;
    Index max_rank = -1;  // <= 0 means "no explicit cap"

    // If true, build U from a frequency-averaged dynamical metric
    //   C_dyn = sum_k Pi(iw_k)^2
    // over the imaginary-axis grid. If false, fall back to Pi(0).
    bool dynamic_basis = true;

    // Hybrid sigma-aware compression:
    //   C = (1-mix) * normalized(C_dyn or C_stat) + mix * normalized(C_sigma)
    //
    // where
    //   C_sigma = sum_{i in QP window} sum_m c_im c_im^T
    // with c_im the full auxiliary-space coupling vector for target QP state i
    // and intermediate state m.
    bool sigma_aware_basis = false;
    double sigma_mix = 0.25;  // in [0,1]

    // Normalize metric blocks before mixing so that the relative weight is
    // controlled by sigma_mix rather than raw matrix magnitude.
    bool normalize_metric_components = true;

    // Print basis construction diagnostics.
    bool print_basis_diagnostics = true;
  };

  struct wc_diagnostic {
    double omega = 0.0;
    double norm_full = 0.0;
    double norm_proj = 0.0;
    double abs_diff = 0.0;
    double rel_diff = 0.0;
  };

  explicit RPA_RI_Reduced(const TCMatrix_gwbse& Mmn) : Mmn_(Mmn) {}

  void configure(Index homo, Index rpamin, Index rpamax) {
    homo_ = homo;
    rpamin_ = rpamin;
    rpamax_ = rpamax;
  }

  void configure_qp_window(Index qpmin, Index qpmax) {
    qpmin_ = qpmin;
    qpmax_ = qpmax;
  }

  void configure_reduced(options opt) { opt_red_ = opt; }

  void setRPAInputEnergies(const Eigen::VectorXd& energies) { energies_ = energies; }
  const Eigen::VectorXd& getRPAInputEnergies() const { return energies_; }

  Eigen::MatrixXd BuildPiImag(double omega) const;
  void BuildReducedBasis();

  Eigen::MatrixXd BuildReducedPiImag(double omega) const;
  Eigen::MatrixXd BuildReducedWcImag(double omega) const;

  Eigen::MatrixXd BuildWcImag(double omega) const;
  Eigen::MatrixXd BuildProjectedReducedWcImag(double omega) const;
  wc_diagnostic CompareWcImag(double omega) const;

  const Eigen::MatrixXd& BasisU() const { return U_; }
  Index rank() const { return U_.cols(); }

 private:
  struct RPAWindow {
    Index lumo;
    Index n_occ;
    Index n_unocc;
  };

  RPAWindow GetWindow() const;
  std::vector<double> ImagFrequencyGrid() const;

  Eigen::MatrixXd BuildStaticOrDynamicCompressionMatrix() const;
  Eigen::MatrixXd BuildSigmaAwareCompressionMatrix() const;
  Eigen::MatrixXd BuildCompressionMatrix() const;

  static Eigen::MatrixXd NormalizeMetric(const Eigen::MatrixXd& M);

  const TCMatrix_gwbse& Mmn_;

  Index homo_ = 0;
  Index rpamin_ = 0;
  Index rpamax_ = 0;

  Index qpmin_ = 0;
  Index qpmax_ = -1;

  Eigen::VectorXd energies_;
  options opt_red_;
  Eigen::MatrixXd U_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_RPA_RI_REDUCED_H