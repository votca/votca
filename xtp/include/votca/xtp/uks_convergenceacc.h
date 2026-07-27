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
#include "votca/xtp/davidsonsolver.h"
#include "votca/xtp/diis.h"
#include "votca/xtp/logger.h"

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

  /// Coupled analogue of BuildSigmaVector: rotates BOTH spin channels
  /// SIMULTANEOUSLY by +-a small step (unlike BuildSigmaVector, which
  /// rotates one channel while holding the other fixed), builds BOTH
  /// new densities together, and calls coupled_fock_builder ONCE per
  /// perturbation to get BOTH new Fock matrices together -- capturing
  /// the genuine alpha-beta coupling (through the shared Coulomb
  /// potential and the XC kernel's cross-spin terms) that
  /// BuildSigmaVector's own, deliberately decoupled treatment
  /// discards. Returns the combined sigma vector (alpha's own block
  /// first, beta's immediately after, matching v's own layout).
  /// Central difference, for the identical reason BuildSigmaVector
  /// itself uses one (see that function's own header comment).
  Eigen::VectorXd BuildCoupledSigmaVector(
      const Eigen::VectorXd& v, const Eigen::MatrixXd& C_alpha,
      Index nocclevels_alpha, const Eigen::MatrixXd& C_beta,
      Index nocclevels_beta, const CoupledFockBuilder& coupled_fock_builder,
      double finite_diff_step = 1e-3) const;

  /// Solves ONE augmented-Hessian eigenvalue problem over the COMBINED
  /// (alpha+beta) rotation space, rather than treating the two spin
  /// channels independently -- see this class's own CoupledFockBuilder
  /// comment for the full motivation. Returns BOTH new MO coefficient
  /// matrices together, since a single, coupled step inherently
  /// produces both simultaneously. predicted_energy_change (output):
  /// the ONE, combined quadratic model's own predicted energy change
  /// for this step -- Iterate's own Fletcher accept/reject logic uses
  /// this directly.
  std::pair<Eigen::MatrixXd, Eigen::MatrixXd> CoupledAugmentedHessianStep(
      const Eigen::MatrixXd& H_AO_alpha, const tools::EigenSystem& MOs_alpha,
      Index nocclevels_alpha, const Eigen::MatrixXd& H_AO_beta,
      const tools::EigenSystem& MOs_beta, Index nocclevels_beta,
      const CoupledFockBuilder& coupled_fock_builder, double trust_radius,
      double& predicted_energy_change) const;

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

  /// Approximate ARH-style direct energy minimization step: builds new
  /// MO coefficients directly via an orbital rotation (occ-virt block
  /// only, since occ-occ/virt-virt rotations leave the density matrix
  /// unchanged), using the orbital gradient (the occ-virt block of the
  /// MO-basis Fock matrix, which vanishes exactly at self-consistency)
  /// and an approximate, diagonal Hessian (orbital energy differences --
  /// the same cheap approximation used as the starting point for most
  /// quasi-Newton SCF methods, e.g. SOSCF's own initial Hessian guess).
  /// The step is trust-radius bounded, and the resulting (only
  /// approximately unitary) rotation is re-orthonormalized via the same
  /// Lowdin approach as OrthogonalizeGuess. See the design discussion
  /// this grew out of for the full reasoning: this is a deliberately
  /// simplified relative of the ARH/TRAH family of direct-minimization
  /// SCF methods (confirmed directly, on this exact system, to resolve
  /// an ADIIS/DIIS convergence failure that ORCA's own auto-TRAH
  /// fallback also had to invoke), not a full implementation of either.
  /// predicted_energy_change (output): the quadratic model's own
  /// predicted E_new - E_old, Sum_ia[g_ia*kappa_ia + 0.5*h_ia*kappa_ia^2]
  /// using the SAME (possibly per-element- and trust-radius-clamped)
  /// kappa actually applied -- needed by Iterate's own Fletcher-style
  /// accept/reject logic (see its own header comment) to compute the
  /// ratio r = E_actual_change / E_predicted_change once the actual,
  /// post-step energy becomes available (on the NEXT Iterate() call,
  /// once the caller has built a new Fock matrix from this step's own
  /// density and passed its energy back in).
  ///
  /// STATUS: superseded by CoupledAugmentedHessianStep below whenever a
  /// coupled Fock-builder callback has been injected (i.e. DFTEngine
  /// has wired up setCoupledFockBuilder) -- kept as the fallback for
  /// any caller that has not done so (e.g.
  /// DFTEngine::RunAtomicDFT_unrestricted), and because
  /// CoupledAugmentedHessianStep itself still uses this same diagonal
  /// approximation as the Davidson preconditioner, exactly as the
  /// TRAH paper's own Sec. II B describes doing.
  Eigen::MatrixXd DirectMinimizationRotation(
      const Eigen::MatrixXd& H_AO, const tools::EigenSystem& MOs,
      Index nocclevels, double& predicted_energy_change) const;

  options opt_alpha_;
  options opt_beta_;

  AOOverlap* S_ = nullptr;
  Logger* log_ = nullptr;
  Eigen::MatrixXd Sminusahalf;

  CoupledFockBuilder coupled_fock_builder_;

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
  double maxerror_ = -1.0;
  Index maxerrorindex_ = 0;
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

  // Fletcher's trust-radius update (Helmich-Paris, J. Chem. Phys. 154,
  // 164104 (2021), Sec. II D -- ORCA's own TRAH-SCF paper, confirmed
  // directly by reading it rather than reconstructed from memory): the
  // trust radius is no longer a fixed constant. direct_min_pending_
  // marks that the MOs/density just returned came from a direct-
  // minimization step (either DirectMinimizationRotation or
  // CoupledAugmentedHessianStep) whose actual effect on the energy has
  // not yet been verified -- checked at the START of the NEXT Iterate()
  // call (see this class's own header comment on
  // DirectMinimizationRotation for why it can only be checked then,
  // not within the same call that took the step).
  bool direct_min_pending_ = false;
  double direct_min_pre_energy_ = 0.0;
  double direct_min_predicted_change_ = 0.0;
  Eigen::MatrixXd direct_min_pre_MOs_alpha_;
  Eigen::MatrixXd direct_min_pre_MOs_beta_;
  Eigen::VectorXd direct_min_pre_MOs_alpha_energies_;
  Eigen::VectorXd direct_min_pre_MOs_beta_energies_;
  double trust_radius_current_ = 0.2;

  // Trust-radius floor tied to BuildSigmaVector's own finite-
  // difference step size (kFiniteDiffStep = 1e-3, defined there): a
  // real run showed trust_radius_current_ shrinking past 1e-6 with
  // an IDENTICAL step/predicted-change/actual-change every single
  // time, since the bisection's own alpha_min=1 floor means the
  // gentlest achievable step cannot shrink any further -- the reject
  // loop could never resolve on its own and only ended when the
  // outer SCF's own iteration budget ran out. A trust radius below
  // the sigma vector's own probing resolution is also asking for
  // precision the underlying finite-difference model was never built
  // to provide. direct_min_floor_hit_ marks that this floor has been
  // reached without an accepted step, at which point Iterate() stops
  // retriggering CoupledAugmentedHessianStep/DirectMinimizationRotation
  // and falls back to plain mixing instead -- reset naturally each time a
  // new UKSConvergenceAcc is constructed (a fresh instance per
  // DFTEngine::EvaluateUKS call, i.e. per outer CDFT lambda trial or
  // per geometry step), not explicitly reset within one instance's
  // own lifetime.
  bool direct_min_floor_hit_ = false;
  // Lowered from 1e-2 to 3e-3 -- still 3x above BuildSigmaVector's own
  // finite-difference resolution (kFiniteDiffStep = 1e-3), so the
  // original justification for having a floor at all still holds, but
  // a real run (He2+, the coupled alpha-beta path) showed
  // actual_dE shrinking monotonically and r steadily climbing back
  // toward the accept range (r=0) right up until it hit the OLD
  // 1e-2 floor and fell back to mixing -- which then settled into a
  // 2-period oscillation rather than converging. Giving that
  // trajectory more room, rather than cutting it off exactly where it
  // was still improving, is a direct, testable next step.
  static constexpr double kMinTrustRadius = 3e-3;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_UKS_CONVERGENCEACC_H