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
#ifndef VOTCA_XTP_CONVERGENCEACC_H
#define VOTCA_XTP_CONVERGENCEACC_H

// VOTCA includes
#include <votca/tools/linalg.h>

// Local VOTCA includes
#include "adiis.h"
#include "aomatrix.h"
#include "diis.h"
#include "logger.h"

namespace votca {
namespace xtp {

/**
 * SCF convergence accelerator for Kohn-Sham iterations.
 *
 * The class stores the Fock and density history needed for linear mixing,
 * Pulay DIIS, and ADIIS, and provides the density construction matching the
 * selected occupation model.
 */
class ConvergenceAcc {
 public:
  /// Occupation model used when constructing density matrices.
  enum KSmode { closed, open, fractional, restricted_open };

  /// User-configurable settings controlling mixing, DIIS, and convergence
  /// thresholds.
  struct options {
    KSmode mode = KSmode::closed;
    bool usediis;
    bool noisy = false;
    Index histlength;
    bool maxout;
    double adiis_start;
    double diis_start;
    double levelshift;
    double levelshiftend;
    // Independent from adiis_start deliberately -- see
    // UKSConvergenceAcc::Iterate's own comment on where this is used
    // for the full reasoning: ORCA's own DampErr is kept fully
    // separate from its own DIISStart, and their own guidance for
    // difficult systems is to make DampErr much SMALLER than default
    // (keeping damping active LONGER), independent of when DIIS
    // itself engages. Reusing adiis_start for both purposes could not
    // represent that independently.
    double mixingend;
    Index numberofelectrons;
    double mixingparameter;
    // Ceiling mixingparameter can adaptively ramp up toward as the SCF
    // struggles, rather than staying fixed at mixingparameter (the
    // BASE/starting value) for the whole run -- matches ORCA's own
    // static-damping design directly (confirmed from a real ORCA log's
    // own resolved SCF settings, not the manual's generic defaults):
    // DampFac (the base, 0.7 by default) and DampMax (the ceiling, 0.98
    // by default) are two separate parameters there, not one fixed
    // value. Notably, 0.98 is exactly the value this session's own
    // hand-tuning independently landed on for a difficult water dimer
    // case -- this adaptive design is a more principled way to obtain
    // that same benefit only when actually needed, rather than paying
    // its cost (slower convergence on easy iterations) for an entire
    // run regardless of whether the system is struggling at all.
    double mixingmax;
    double Econverged;
    double error_converged;
    Index number_alpha_electrons = 0;
    Index number_beta_electrons = 0;
    // Maximum iterations for the Davidson eigensolver used by
    // CoupledAugmentedHessianStep's own direct-minimization fallback --
    // NOT CDFT-specific, since that fallback can engage for any
    // sufficiently difficult UKS SCF (see this field's own XML help
    // text, dftpackage.xml, for the real case that motivated exposing
    // this at all: a strong CDFT constraint over a large fragment left
    // the solver still short of its own convergence tolerance at the
    // previous, hardcoded default of 50).
    Index davidson_max_iter = 50;
    // Energy rise (Hartree) above the lowest energy reached so far at which
    // the SCF discards its extrapolation history and restarts from the
    // lowest-energy density. 0 disables.
    double energy_reset = 1.0;
  };

  /// Spin-resolved density matrices returned for open-shell SCF updates.
  struct SpinDensity {
    Eigen::MatrixXd alpha;
    Eigen::MatrixXd beta;

    /// Return the total density P = P^alpha + P^beta.
    Eigen::MatrixXd total() const { return alpha + beta; }

    /// Return the spin density P^alpha - P^beta.
    Eigen::MatrixXd spin() const { return alpha - beta; }
  };

  /// Store SCF acceleration settings and derive the number of occupied levels
  /// for the selected KS mode.
  void Configure(const ConvergenceAcc::options& opt) {
    opt_ = opt;
    if (opt_.mode == KSmode::closed) {
      nocclevels_ = opt_.numberofelectrons / 2;
    } else if (opt_.mode == KSmode::open) {
      nocclevels_ = opt_.numberofelectrons;
    } else if (opt_.mode == KSmode::fractional) {
      nocclevels_ = 0;
    } else if (opt_.mode == KSmode::restricted_open) {
      nocclevels_ =
          std::max(opt_.number_alpha_electrons, opt_.number_beta_electrons);
    }
    diis_.setHistLength(opt_.histlength);
    StartNewSCF();
  }

  // Forget the lowest-energy point used by the energy-rise reset. Called
  // at the start of every SCF: successive SCFs on one accelerator (e.g.
  // the stages of DFT-in-DFT embedding) need not share an energy scale.
  // The extrapolation history itself is left as it was.
  void StartNewSCF() {
    have_best_ = false;
    energy_resets_ = 0;
    iterations_since_reset_ = kResetCooldown;
  }
  /// Attach the logger used for convergence diagnostics.
  void setLogger(Logger* log) { log_ = log; }

  /// Print the active convergence-acceleration settings to the logger.
  void PrintConfigOptions() const;

  /// Check whether both the total-energy change and DIIS error are below their
  /// thresholds.
  bool isConverged() const {
    if (totE_.size() < 2) {
      return false;
    } else {
      return std::abs(getDeltaE()) < opt_.Econverged &&
             getDIIsError() < opt_.error_converged;
    }
  }

  /// Return the total-energy change between the two most recent SCF iterations.
  double getDeltaE() const {
    if (totE_.size() < 2) {
      return 0;
    } else {
      return totE_.back() - totE_[totE_.size() - 2];
    }
  }
  /// Precompute overlap-dependent quantities used when solving the Fock matrix.
  // Builds X = S^-1/2 for the orthogonal basis the Fock matrix is
  // diagonalized in. Eigenvalues of S below etol are dropped from X; the
  // directions they span are pushed to the top of every Fock spectrum
  // (kRemovedShift) so that they can never be occupied or appear among the
  // low virtuals. Their MO coefficient vectors are zero.
  void setOverlap(AOOverlap& S, double etol);

  /// Return the DIIS commutator norm from the latest iteration.
  double getDIIsError() const { return diiserror_; }

  /// Report whether plain density mixing is currently used instead of
  /// extrapolation.
  bool getUseMixing() const { return usedmixing_; }

  // Consistency check of the extrapolation history: every stored
  // (Fock, density) pair must still produce the error matrix DIIS holds at
  // the same position. Used by the unit tests.
  bool HistoryIsAligned() const;

  /// Advance the SCF accelerator by one step and return the updated density
  /// matrix.
  Eigen::MatrixXd Iterate(const Eigen::MatrixXd& dmat, Eigen::MatrixXd& H,
                          tools::EigenSystem& MOs, double totE);
  /// Solve the generalized eigenvalue problem for the current Fock matrix.
  tools::EigenSystem SolveFockmatrix(const Eigen::MatrixXd& H) const;
  /// Apply a virtual-space level shift in the molecular-orbital basis.
  void Levelshift(Eigen::MatrixXd& H, const Eigen::MatrixXd& MOs_old) const;

  /// Build the density matrix corresponding to the configured KS occupation
  /// model.
  Eigen::MatrixXd DensityMatrix(const tools::EigenSystem& MOs) const;

  /// Build separate alpha and beta density matrices for spin-resolved SCF
  /// modes.
  SpinDensity DensityMatrixSpinResolved(const tools::EigenSystem& MOs) const;

 private:
  options opt_;

  /// Construct a closed-shell ground-state density matrix from occupied
  /// orbitals.
  Eigen::MatrixXd DensityMatrixGroundState(const Eigen::MatrixXd& MOs) const;
  /// Construct a fully occupied unrestricted density matrix from the supplied
  /// orbitals.
  Eigen::MatrixXd DensityMatrixGroundState_unres(
      const Eigen::MatrixXd& MOs) const;
  /// Construct a fractional-occupation density matrix from orbital occupations.
  Eigen::MatrixXd DensityMatrixGroundState_frac(
      const tools::EigenSystem& MOs) const;

  /// Construct alpha and beta densities for a restricted open-shell
  /// determinant.
  SpinDensity DensityMatrixGroundState_restricted_open(
      const Eigen::MatrixXd& MOs) const;

  bool usedmixing_ = true;
  double diiserror_ = std::numeric_limits<double>::max();
  Logger* log_;
  const AOOverlap* S_;

  Eigen::MatrixXd Sminusahalf;
  // Projector onto the directions removed from Sminusahalf, in the
  // orthogonal coordinates of SolveFockmatrix; empty if none were removed.
  Eigen::MatrixXd removed_projector_;
  static constexpr double kRemovedShift = 1e3;  // Hartree

  // Extrapolation history, all three aligned entry by entry: Fock matrix,
  // the density it was built from, and that pair's DIIS error. DIIS keeps
  // its own copy of the error matrices and is trimmed at the same index.
  std::vector<Eigen::MatrixXd> mathist_;
  std::vector<Eigen::MatrixXd> dmatHist_;
  std::vector<double> errhist_;
  // Every energy, in order. Not trimmed with the history: DeltaE and the
  // energy tests compare consecutive iterations.
  std::vector<double> totE_;

  // Lowest-energy point so far, for the energy-rise reset.
  bool have_best_ = false;
  double best_energy_ = 0.0;
  Eigen::MatrixXd best_dmat_;
  Eigen::MatrixXd best_H_;
  Index energy_resets_ = 0;
  // A reset may not fire again until the history has rebuilt, and only a
  // few times per SCF: returning to the same point over and over must be
  // impossible.
  static constexpr Index kResetCooldown = 3;
  static constexpr Index kMaxEnergyResets = 5;
  Index iterations_since_reset_ = kResetCooldown;

  // ADIIS instead of DIIS once the energy rises by more than this
  // (Hartree) between iterations, even below DIIS_start.
  static constexpr double kEnergyRiseForADIIS = 1e-4;

  Index nocclevels_;
  ADIIS adiis_;
  DIIS diis_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CONVERGENCEACC_H
