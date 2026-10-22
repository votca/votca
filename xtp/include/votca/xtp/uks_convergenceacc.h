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

#ifndef VOTCA_XTP_UKS_CONVERGENCEACC_H
#define VOTCA_XTP_UKS_CONVERGENCEACC_H

// Standard includes
#include <functional>
#include <utility>

// VOTCA includes
#include <votca/tools/linalg.h>

#include "votca/xtp/adiis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/convergenceacc.h"
#include "votca/xtp/diis.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/trah.h"

namespace votca {
namespace xtp {

class UKSConvergenceAcc {
 public:
  using options = ConvergenceAcc::options;
  using KSmode = ConvergenceAcc::KSmode;

  struct SpinDensity {
    Eigen::MatrixXd alpha;
    Eigen::MatrixXd beta;

    Eigen::MatrixXd total() const { return alpha + beta; }
  };

  struct SpinFock {
    Eigen::MatrixXd alpha;
    Eigen::MatrixXd beta;
  };

  /// Given BOTH new spin densities together, returns BOTH new AO Fock
  /// matrices together (H0 + Coulomb + exchange + XC, for each
  /// channel). Closer to VOTCA's own native Fock-build structure:
  /// Vxc_Potential::IntegrateVXCSpin already computes vxc_alpha and
  /// vxc_beta together from both densities in one call. Used by
  /// CoupledAugmentedHessianStep/BuildCoupledSigmaVector to capture
  /// the genuine alpha-beta coupling (through the shared Coulomb
  /// potential and the XC kernel's cross-spin terms) that treating the
  /// two channels independently would discard -- see the conversation
  /// this grew out of: a direct, controlled comparison against ORCA on
  /// an identical geometry showed ORCA's own, fully-coupled TRAH
  /// converging where an earlier, decoupled implementation here did
  /// not. Injected from DFTEngine (setCoupledFockBuilder below) since
  /// the actual integral/XC machinery (CalcERIs/CalcERIs_EXX,
  /// Vxc_Potential::IntegrateVXCSpin) lives there, not in this class --
  /// matching the same pointer-injection pattern already used for
  /// S_/log_.
  using CoupledFockBuilder =
      std::function<SpinFock(const Eigen::MatrixXd& /*Dalpha_new*/,
                             const Eigen::MatrixXd& /*Dbeta_new*/)>;

  void Configure(const options& opt_alpha, const options& opt_beta);
  void setLogger(Logger* log);
  void setOverlap(AOOverlap& S, double etol);
  void setCoupledFockBuilder(const CoupledFockBuilder& builder) {
    coupled_fock_builder_ = builder;
  }

  /// Analytic UKS energy Hessian d^2E/dx dx (both channels, layout of
  /// trah.h) at the orbitals (MOs_alpha, MOs_beta) with Fock matrices
  /// (F_alpha, F_beta), applied to the columns of X. If set, it replaces
  /// the finite-difference products of BuildCoupledSigmaVector; if it
  /// throws (no second derivatives in libxc), those are used instead.
  using OrbitalHessian = std::function<Eigen::MatrixXd(
      const tools::EigenSystem& MOs_alpha, const tools::EigenSystem& MOs_beta,
      const Eigen::MatrixXd& F_alpha, const Eigen::MatrixXd& F_beta,
      const Eigen::MatrixXd& X)>;
  void setOrbitalHessian(const OrbitalHessian& hessian) {
    orbital_hessian_ = hessian;
  }

  SpinDensity DensityMatrix(const tools::EigenSystem& MOs_alpha,
                            const tools::EigenSystem& MOs_beta) const;

  SpinDensity Iterate(const SpinDensity& dmat, SpinFock& H,
                      tools::EigenSystem& MOs_alpha,
                      tools::EigenSystem& MOs_beta, double totE);

  tools::EigenSystem SolveFockmatrix(const Eigen::MatrixXd& H) const;

  /// One column of the Davidson trial-vector matrix V, containing
  /// [v0; v_ov] (the scalar augmented component and the flattened
  /// occ-virt rotation-space vector), mapped back to a full nao x nao
  /// antisymmetric matrix (occ-virt/virt-occ blocks only, matching
  /// DirectMinimizationRotation's own convention) -- needed to
  /// actually apply a trial rotation direction when building the
  /// finite-difference sigma vector. Public so it can be exercised
  /// directly by unit tests (see test_uks_convergenceacc.cc), given how
  /// cheaply this pure, self-contained math can be checked in isolation
  /// compared to the SCF machinery that otherwise relies on it.
  Eigen::MatrixXd UnflattenRotation(const Eigen::VectorXd& v_ov, Index nao,
                                    Index nocclevels) const;

  /// Coupled analogue of UnflattenRotation: splits ONE combined vector
  /// v (size n_ov_alpha + n_ov_beta -- alpha's own occ-virt block
  /// first, beta's own immediately after) into TWO separate
  /// antisymmetric rotation matrices, kappa_alpha (nao_alpha x
  /// nao_alpha) and kappa_beta (nao_beta x nao_beta). Public for the
  /// same reason as UnflattenRotation itself.
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> UnflattenCoupledRotation(
      const Eigen::VectorXd& v, Index nao_alpha, Index nocclevels_alpha,
      Index nao_beta, Index nocclevels_beta) const;

  /// Orbital rotations C' = C (1 + kappa), with kappa_ia = x_ia (i
  /// occupied, a virtual) and kappa_ai = -x_ia, followed by Lowdin
  /// orthonormalisation. Occupied orbital i changes by -x_ia phi_a, so for
  /// one spin channel the energy gradient is dE/dx_ia = -2 F_ia, with F the
  /// Fock matrix in the current MO basis. All gradients, Hessian products
  /// and predicted energy changes of the second-order fallback use this
  /// scale.
  static constexpr double EnergyGradientScale = -2.0;

  /// Product of the energy Hessian d^2E/dx dx (both spin channels, with
  /// the alpha-beta coupling through J and the XC kernel) with v: a
  /// central finite difference of EnergyGradientScale * F_ia along v, with
  /// BOTH channels rotated together and one coupled_fock_builder call per
  /// displacement. Layout as UnflattenCoupledRotation: alpha block first.
  Eigen::VectorXd BuildCoupledSigmaVector(
      const Eigen::VectorXd& v, const Eigen::MatrixXd& C_alpha,
      Index nocclevels_alpha, const Eigen::MatrixXd& C_beta,
      Index nocclevels_beta, const CoupledFockBuilder& coupled_fock_builder,
      double finite_diff_step = 1e-3) const;

  /// TRAH step (trah::Step) over the combined (alpha+beta) rotation space,
  /// with Hessian products from the analytic OrbitalHessian if set (else
  /// BuildCoupledSigmaVector with coupled_fock_builder) and the orbitals
  /// rotated by exp(kappa). Returns both new MO coefficient matrices and, in
  /// predicted_energy_change, g.x + 1/2 x.H.x of the step taken (used by
  /// Iterate's Fletcher accept/reject check).
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> CoupledAugmentedHessianStep(
      const Eigen::MatrixXd& H_AO_alpha, const tools::EigenSystem& MOs_alpha,
      Index nocclevels_alpha, const Eigen::MatrixXd& H_AO_beta,
      const tools::EigenSystem& MOs_beta, Index nocclevels_beta,
      const CoupledFockBuilder& coupled_fock_builder, double trust_radius,
      double& predicted_energy_change) const;

  /// Diagonal-Hessian direct-minimisation step for one spin channel,
  /// used when no coupled Fock builder is set (e.g. the atomic UKS guess):
  /// kappa_ia = -g_ia / h_ia with g_ia = -2 F_ia and h_ia = 2 (eps_a -
  /// eps_i) (gaps below 1e-3 Ha raised to 1e-3), each element clamped to
  /// 0.1 and the whole step to the current trust radius (Frobenius norm of
  /// the antisymmetric kappa). predicted_energy_change: the model's
  /// Sum_ia [g_ia kappa_ia + 1/2 h_ia kappa_ia^2] for the step taken.
  Eigen::MatrixXd DirectMinimizationRotation(
      const Eigen::MatrixXd& H_AO, const tools::EigenSystem& MOs,
      Index nocclevels, double& predicted_energy_change) const;

  bool isConverged() const;
  double getDIIsError() const { return diiserror_; }
  double getDeltaE() const;
  bool getUseMixing() const { return usedmixing_; }

 private:
  Eigen::MatrixXd DensityMatrixGroundState_unres(const Eigen::MatrixXd& MOs,
                                                 Index nocclevels) const;

  void Levelshift(Eigen::MatrixXd& H, const Eigen::MatrixXd& MOs_old,
                  const options& opt, Index nocclevels) const;

  Eigen::MatrixXd BuildErrorMatrix(const Eigen::MatrixXd& dmat,
                                   const Eigen::MatrixXd& H) const;

  double CombinedError(const Eigen::MatrixXd& err_alpha,
                       const Eigen::MatrixXd& err_beta) const;

  bool TrailingAverageStalled() const;

  /// One second-order step from the current orbitals (or their density
  /// first, if dmat is not it); see Iterate.
  SpinDensity SecondOrderStep(const SpinDensity& dmat, const SpinFock& H,
                              tools::EigenSystem& MOs_alpha,
                              tools::EigenSystem& MOs_beta, double totE);

  options opt_alpha_;
  options opt_beta_;

  AOOverlap* S_ = nullptr;
  Logger* log_ = nullptr;
  Eigen::MatrixXd Sminusahalf;

  CoupledFockBuilder coupled_fock_builder_;
  OrbitalHessian orbital_hessian_;
  mutable bool analytic_hessian_failed_ = false;

  Index nocclevels_alpha_ = 0;
  Index nocclevels_beta_ = 0;

  std::vector<Eigen::MatrixXd> mathist_alpha_;
  std::vector<Eigen::MatrixXd> mathist_beta_;
  std::vector<Eigen::MatrixXd> dmatHist_alpha_;
  std::vector<Eigen::MatrixXd> dmatHist_beta_;
  std::vector<double> totE_;

  DIIS diis_;
  ADIIS adiis_;

  double diiserror_ = 1.0;
  // DIIS error of each history entry, aligned with mathist_alpha_.
  std::vector<double> errhist_;
  // See ConvergenceAcc.
  static constexpr double kEnergyRiseForADIIS = 1e-4;
  Eigen::MatrixXd removed_projector_;
  static constexpr double kRemovedShift = 1e3;
  bool usedmixing_ = true;

  // Direct-minimization fallback bookkeeping. consecutive_adiis_failures_
  // tracks how many (A)DIIS attempts in a row have failed -- switching
  // to DirectMinimizationRotation/CoupledAugmentedHessianStep after
  // kMaxConsecutiveADIISFailures deliberately mirrors ORCA's own
  // auto-TRAH trigger (switching away from DIIS-family methods after
  // they visibly struggle), rather than falling back to plain mixing
  // indefinitely the way this class already did before this addition.
  Index consecutive_adiis_failures_ = 0;
  static constexpr Index kMaxConsecutiveADIISFailures = 5;

  // Trailing-average trigger for direct-minimization, ADDITIONAL to
  // (not replacing) the consecutive-failures count above -- matches
  // ORCA's own AutoTRAH design (confirmed directly from a real ORCA
  // log's own resolved SCF settings: "Auto Start start iteration 50",
  // "Auto Start num. interpolation iter. 10", "Auto Start mean grad.
  // ratio tolernc. 1.125"). Rather than only reacting to N CONSECUTIVE
  // ADIIS failures (which says nothing about whether progress is
  // merely slow but real, or genuinely stalled), this tracks whether
  // diiserror_ itself is failing to shrink fast enough ON AVERAGE over
  // a trailing window, regardless of whether any individual ADIIS
  // attempt in that window happened to "succeed" by its own tail-
  // coefficient criterion. diiserror_history_ holds the trailing
  // window (capped at kTrailingWindowSize entries, oldest dropped as
  // new ones are added); the trigger itself only becomes eligible
  // after kAutoStartIteration total iterations, matching ORCA's own
  // deliberate delay before considering this criterion at all.
  std::vector<double> diiserror_history_;
  Index total_iteration_count_ = 0;
  static constexpr Index kTrailingWindowSize = 10;
  static constexpr Index kAutoStartIteration = 50;
  static constexpr double kMeanRatioTolerance = 1.125;

  // A second-order step has been returned and is checked (Fletcher,
  // trust_) at the start of the next Iterate call, once its energy is known;
  // the pre-step state is kept to undo a rejected step.
  bool direct_min_pending_ = false;
  double direct_min_pre_energy_ = 0.0;
  double direct_min_predicted_change_ = 0.0;
  Eigen::MatrixXd direct_min_pre_MOs_alpha_;
  Eigen::MatrixXd direct_min_pre_MOs_beta_;
  Eigen::VectorXd direct_min_pre_MOs_alpha_energies_;
  Eigen::VectorXd direct_min_pre_MOs_beta_energies_;
  SpinFock direct_min_pre_H_;
  double direct_min_linear_ = 0.0;
  double direct_min_step_norm_ = 0.0;
  mutable double last_step_linear_ = 0.0;
  mutable double last_step_norm_ = 0.0;
  // subspace of the last step, reused when a rejected step is repeated
  mutable trah::Subspace trah_subspace_;
  bool last_rejected_ = false;
  trah::TrustRadius trust_{0.5};
  trah::StepOptions trah_options_;
  // second-order (TRAH) steps in use for this SCF
  bool trah_active_ = false;

  // Trust-radius floor: once a rejected step shrinks the trust radius
  // below it, direct minimisation is given up for this SCF and Iterate
  // falls back to mixing (direct_min_floor_hit_). Reset with each new
  // UKSConvergenceAcc (one per DFTEngine::EvaluateUKS call).
  bool direct_min_floor_hit_ = false;
  static constexpr double kMinTrustRadius = 3e-3;
  // Largest element of |D(MOs) - D| for which the input density is
  // treated as the density of the current orbitals (see Iterate).
  static constexpr double kDensityMismatchTolerance = 1e-8;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_UKS_CONVERGENCEACC_H