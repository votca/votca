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

#pragma once
#ifndef VOTCA_XTP_GW_H
#define VOTCA_XTP_GW_H

#include <unordered_set>

// Local VOTCA includes
#include "logger.h"
#include "orbitals.h"
#include "qp_solver_utils.h"
#include "rpa.h"
#include "sigma_base.h"
#include "threecenter.h"

namespace votca {
namespace xtp {

class GW {

  using EvalStage = qp_solver::EvalStage;
  using QPStats = qp_solver::Stats;
  using QPRootCandidate = qp_solver::RootCandidate;
  using QPWindowDiagnostics = qp_solver::WindowDiagnostics;

 public:
  GW(Logger& log, TCMatrix_gwbse& Mmn, const Eigen::MatrixXd& vxc,
     const Eigen::VectorXd& dft_energies)
      : log_(log),
        Mmn_(Mmn),
        vxc_(vxc),
        dft_energies_(dft_energies),
        rpa_(log, Mmn) {};

  struct options {
    Index homo;
    Index qpmin;
    Index qpmax;
    Index rpamin;
    Index rpamax;
    double eta;
    double g_sc_limit;
    Index g_sc_max_iterations;
    double gw_sc_limit;
    Index gw_sc_max_iterations;
    double shift = 0;
    double ScaHFX = 0.0;
    std::string sigma_integration;
    Index reset_3c;  // how often the 3c integrals in iterate should be
                     // rebuilt
    std::string qp_solver;
    double qp_solver_alpha = 0.75;
    // Legacy aliases from the original grid-search implementation.
    // They used to define both the total QP window and the scan resolution at
    // the same time. The new implementation keeps them only so that old XML
    // files and existing tests still work; configuration code maps them onto
    // the decoupled controls below.
    Index qp_grid_steps = 0;
    double qp_grid_spacing = 0.0;

    // Decoupled QP search controls:
    //  - qp_full_window_half_width defines the outer full search interval
    //  - qp_dense_spacing defines the dense scan used for robust sign changes
    //  - qp_adaptive_shell_width / count define the outward shell search
    // All values are normalized in configure(); negative values mean "unset"
    // and trigger either legacy mapping or robust defaults.
    double qp_full_window_half_width = -1.0;  // Ha; <=0 means "unset"
    double qp_dense_spacing = -1.0;           // Ha; <=0 means "unset"
    double qp_adaptive_shell_width = -1.0;    // Ha; <=0 means "unset"
    Index qp_adaptive_shell_count = 0;        // 0 means "use shell width"

    Index gw_mixing_order;          // mixing order
    double gw_mixing_alpha;         // mixing alpha, also linear mixing
    std::string quadrature_scheme;  // Kind of Gaussian-quadrature scheme to use
    Index order;   // only needed for complex integration sigma CDA
    double alpha;  // smooth tail in complex integration sigma CDA
    bool qp_restrict_search = true;
    double qp_zero_margin = 1e-6;
    double qp_virtual_min_energy = -0.1;
    std::string qp_root_finder = "bisection";
    std::string qp_grid_search_mode = "adaptive_with_dense_fallback";

    // QSGW options
    // When true, run quasiparticle self-consistent GW instead of evGW.
    // The three-centre integrals are rotated at each iteration so that W and
    // Sigma are evaluated in the basis of the current QP wavefunctions.
    bool do_qsgw = false;
    Index qsgw_max_iterations = 20;
    double qsgw_sc_limit = 1e-5;  // Ha; convergence threshold on QP energies

    // Maximum allowed perturbative QP correction for virtual states to be
    // included in the QSGW self-consistency loop. Virtual states with
    // |e_QP - e_DFT| > qsgw_max_virt_correction are excluded from the QSGW
    // window and kept at their perturbative QP energies (DFT-MO wavefunctions).
    // This prevents pathological high-energy basis-set artefact states from
    // destabilising the QSGW rotation. Occupied states are never excluded.
    // Default: 0.5 Ha (13.6 eV, 1 Rydberg) -- corrections above this are
    // unlikely to be physically meaningful for typical molecular systems.
    // Set to a large value (e.g. 1e10) to disable the threshold entirely.
    double qsgw_max_virt_correction = 0.5;  // Ha
  };

  void configure(const options& opt);

  Eigen::VectorXd getGWAResults() const;
  // Calculates the diagonal elements up to self consistency
  void CalculateGWPerturbation();

  // Calculated offdiagonal elements as well
  void CalculateHQP();

  /**
   * \brief Run quasiparticle self-consistent GW (QSGW).
   *
   * At each iteration:
   *   1. Compute the full off-diagonal self-energy matrix tilde_Sigma
   *      (symmetrised static approximation) in the current QP basis.
   *   2. Diagonalise H_QSGW = diag(e_DFT - v_xc) + tilde_Sigma to obtain
   *      new QP energies and a rotation matrix U.
   *   3. Rotate the three-centre integrals Mmn via U so that W and Sigma
   *      are computed from the updated QP wavefunctions next iteration.
   * Converges when max|e_QP^{k+1} - e_QP^k| < qsgw_sc_limit.
   *
   * The accumulated rotation U (DFT MOs -> converged QP wavefunctions) is
   * stored in qsgw_rotation_ and returned via getQSGWRotation().
   * QPdiag energies and eigenvectors are set to the converged QSGW values.
   */
  void CalculateQSGW();

  /// Return the accumulated QSGW rotation matrix U (DFT MOs -> QP
  /// wavefunctions)
  const Eigen::MatrixXd& getQSGWRotation() const { return qsgw_rotation_; }

  /// Return the seed (G0W0/evGW) energies that QSGW started from
  const Eigen::VectorXd& getQSGWSeedEnergies() const {
    return qsgw_seed_energies_;
  }

  /**
   * \brief Print a two-column comparison of seed (G0W0 or evGW) vs converged
   *        QSGW quasiparticle energies.
   *
   * @param seed_label    Label for the seed column, e.g. "G0W0" or "evGW".
   * @param seed_energies QP energies from the seed calculation (Ha).
   * @param qsgw_energies Converged QSGW eigenvalues (Ha).
   */
  void PrintQSGW_Energies(const std::string& seed_label,
                          const Eigen::VectorXd& seed_energies,
                          const Eigen::VectorXd& qsgw_energies) const;

  /**
   * \brief Print the dominant DFT-KS orbital contributions to each converged
   *        QSGW quasiparticle state.
   *
   * For each QP state n, prints the KS orbitals m with |U_{mn}|^2 > threshold,
   * where U = qsgw_rotation_ (columns = QP states in DFT-KS basis).
   *
   * @param threshold  Minimum weight to print (default 0.01 = 1%)
   */
  void PrintQSGW_Composition(double threshold = 0.01) const;

  Eigen::MatrixXd getHQP() const;

  // Diagonalize QP particle Hamiltonian
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagonalizeQPHamiltonian()
      const;

  void PlotSigma(std::string filename, Index steps, double spacing,
                 std::string states) const;

  Eigen::VectorXd RPAInputEnergies() const {
    return rpa_.getRPAInputEnergies();
  }

 private:
  Index qptotal_;

  Eigen::MatrixXd Sigma_x_;
  Eigen::MatrixXd Sigma_c_;
  Eigen::MatrixXd qsgw_rotation_;       // accumulated U: DFT MOs -> QSGW QP
                                        // wavefunctions
  Eigen::VectorXd qsgw_seed_energies_;   // evGW/G0W0 energies used as QSGW seed
  // Merged QSGW+seed energies for the full QP window. Set by CalculateQSGW()
  // when the virtual window is trimmed. getGWAResults() returns this when
  // non-empty, bypassing the sigma-matrix recomputation which would be wrong
  // for the excluded levels.
  Eigen::VectorXd qsgw_final_energies_;

  options opt_;

  std::unique_ptr<Sigma_base> sigma_ = nullptr;
  Logger& log_;
  TCMatrix_gwbse& Mmn_;
  const Eigen::MatrixXd& vxc_;
  const Eigen::VectorXd& dft_energies_;

  Index gw_sc_iteration_;

  RPA rpa_;
  // small class which calculates f(w) with and df/dw(w)
  // f=Sigma_c(w)+offset-w
  // offset= e_dft+Sigma_x-Vxc
  class QPFunc {
   public:
    QPFunc(Index gw_level, const Sigma_base& sigma, double offset)
        : gw_level_(gw_level), offset_(offset), sigma_c_func_(sigma) {}

    std::pair<double, double> operator()(double frequency) const {
      std::pair<double, double> result;
      result.first = value(frequency, EvalStage::Other);
      result.second = deriv(frequency);
      return result;
    }

    double sigma(double frequency, EvalStage stage = EvalStage::Other) const {
      const std::uint64_t key = FrequencyKey(frequency);

      auto insert_result = seen_frequencies_.insert(key);
      if (!insert_result.second) {
        ++stats_.sigma_repeat_calls;
      } else {
        ++stats_.sigma_unique_frequencies;
      }

      CountSigmaStage(stage);
      return sigma_c_func_.CalcCorrelationDiagElement(gw_level_, frequency);
    }

    double value(double frequency, EvalStage stage = EvalStage::Other) const {
      return sigma(frequency, stage) + offset_ - frequency;
    }

    double deriv(double frequency) const {
      ++stats_.deriv_calls;
      return sigma_c_func_.CalcCorrelationDiagElementDerivative(gw_level_,
                                                                frequency) -
             1.0;
    }

    const QPStats& GetStats() const { return stats_; }

   private:
    static std::uint64_t FrequencyKey(double x) {
      std::uint64_t key = 0;
      static_assert(sizeof(double) == sizeof(std::uint64_t),
                    "Unexpected double size");
      std::memcpy(&key, &x, sizeof(double));
      return key;
    }

    void CountSigmaStage(EvalStage stage) const {
      switch (stage) {
        case EvalStage::Scan:
          ++stats_.sigma_scan_calls;
          break;
        case EvalStage::Refine:
          ++stats_.sigma_refine_calls;
          break;
        case EvalStage::Derivative:
          ++stats_.sigma_derivative_calls;
          break;
        case EvalStage::Other:
        default:
          ++stats_.sigma_other_calls;
          break;
      }
    }

    Index gw_level_;
    double offset_;
    const Sigma_base& sigma_c_func_;

    mutable std::unordered_set<std::uint64_t> seen_frequencies_;
    mutable QPStats stats_;
  };

  double CalcHomoLumoShift(Eigen::VectorXd frequencies) const;
  Eigen::VectorXd ScissorShift_DFTlevel(
      const Eigen::VectorXd& dft_energies) const;
  void PrintQP_Energies(const Eigen::VectorXd& qp_diag_energies) const;
  void PrintGWA_Energies() const;

  Eigen::VectorXd SolveQP(const Eigen::VectorXd& frequencies) const;
  boost::optional<double> SolveQP_Grid(double intercept0, double frequency0,
                                       Index gw_level,
                                       QPStats* stats = nullptr) const;

  boost::optional<double> SolveQP_Grid_Windowed(
      double intercept0, double frequency0, Index gw_level, double left_limit,
      double right_limit, bool allow_rejected_return = true,
      QPStats* stats = nullptr) const;

  boost::optional<double> SolveQP_Grid_Windowed_Adaptive(
      double intercept0, double frequency0, Index gw_level, double left_limit,
      double right_limit, bool allow_rejected_return = true,
      QPStats* stats = nullptr) const;

  boost::optional<double> SolveQP_Grid_Windowed_Dense(
      double intercept0, double frequency0, Index gw_level, double left_limit,
      double right_limit, bool allow_rejected_return = true,
      QPStats* stats = nullptr) const;

  boost::optional<double> SolveQP_FixedPoint(double intercept0,
                                             double frequency0, Index gw_level,
                                             QPStats* stats = nullptr) const;
  boost::optional<double> SolveQP_Linearisation(double intercept0,
                                                double frequency0,
                                                Index gw_level,
                                                QPStats* stats = nullptr) const;
  bool Converged(const Eigen::VectorXd& e1, const Eigen::VectorXd& e2,
                 double epsilon) const;

  boost::optional<QPRootCandidate> RefineQPInterval(
      double lowerbound, double f_lowerbound, double upperbound,
      double f_upperbound, const QPFunc& f, double reference) const;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GW_H
