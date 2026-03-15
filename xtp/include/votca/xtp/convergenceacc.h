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

class ConvergenceAcc {
 public:
enum KSmode { closed, open, fractional, restricted_open };

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
    Index numberofelectrons;
    double mixingparameter;
    double Econverged;
    double error_converged;
    Index number_alpha_electrons = 0;
    Index number_beta_electrons = 0;
  };

struct SpinDensity {
  Eigen::MatrixXd alpha;
  Eigen::MatrixXd beta;

  Eigen::MatrixXd total() const { return alpha + beta; }
  Eigen::MatrixXd spin() const { return alpha - beta; }
};

  /// Store SCF acceleration settings and derive the number of occupied levels for the selected KS mode.
  void Configure(const ConvergenceAcc::options& opt) {
    opt_ = opt;
    if (opt_.mode == KSmode::closed) {
      nocclevels_ = opt_.numberofelectrons / 2;
    } else if (opt_.mode == KSmode::open) {
      nocclevels_ = opt_.numberofelectrons;
    } else if (opt_.mode == KSmode::fractional) {
      nocclevels_ = 0;
    } else if (opt_.mode == KSmode::restricted_open) {
      nocclevels_ = std::max(opt_.number_alpha_electrons,
                         opt_.number_beta_electrons);
    }
    diis_.setHistLength(opt_.histlength);
  }
  /// Attach the logger used for convergence diagnostics.
  void setLogger(Logger* log) { log_ = log; }

  /// Print the active convergence-acceleration settings to the logger.
  void PrintConfigOptions() const;

  /// Check whether both the total-energy change and DIIS error are below their thresholds.
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
  void setOverlap(AOOverlap& S, double etol);

  /// Return the DIIS commutator norm from the latest iteration.
  double getDIIsError() const { return diiserror_; }

  /// Report whether plain density mixing is currently used instead of extrapolation.
  bool getUseMixing() const { return usedmixing_; }

  /// Advance the SCF accelerator by one step and return the updated density matrix.
  Eigen::MatrixXd Iterate(const Eigen::MatrixXd& dmat, Eigen::MatrixXd& H,
                          tools::EigenSystem& MOs, double totE);
  /// Solve the generalized eigenvalue problem for the current Fock matrix.
  tools::EigenSystem SolveFockmatrix(const Eigen::MatrixXd& H) const;
  /// Apply a virtual-space level shift in the molecular-orbital basis.
  void Levelshift(Eigen::MatrixXd& H, const Eigen::MatrixXd& MOs_old) const;

  /// Build the density matrix corresponding to the configured KS occupation model.
  Eigen::MatrixXd DensityMatrix(const tools::EigenSystem& MOs) const;

  /// Build separate alpha and beta density matrices for spin-resolved SCF modes.
  SpinDensity DensityMatrixSpinResolved(const tools::EigenSystem& MOs) const;

 private:
  options opt_;

  /// Construct a closed-shell ground-state density matrix from occupied orbitals.
  Eigen::MatrixXd DensityMatrixGroundState(const Eigen::MatrixXd& MOs) const;
  /// Construct a fully occupied unrestricted density matrix from the supplied orbitals.
  Eigen::MatrixXd DensityMatrixGroundState_unres(
      const Eigen::MatrixXd& MOs) const;
  /// Construct a fractional-occupation density matrix from orbital occupations.
  Eigen::MatrixXd DensityMatrixGroundState_frac(
      const tools::EigenSystem& MOs) const;

  /// Construct alpha and beta densities for a restricted open-shell determinant.
  SpinDensity DensityMatrixGroundState_restricted_open(
    const Eigen::MatrixXd& MOs) const;

  bool usedmixing_ = true;
  double diiserror_ = std::numeric_limits<double>::max();
  Logger* log_;
  const AOOverlap* S_;

  Eigen::MatrixXd Sminusahalf;
  std::vector<Eigen::MatrixXd> mathist_;
  std::vector<Eigen::MatrixXd> dmatHist_;
  std::vector<double> totE_;

  Index nocclevels_;
  Index maxerrorindex_ = 0;
  double maxerror_ = 0.0;
  ADIIS adiis_;
  DIIS diis_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CONVERGENCEACC_H
