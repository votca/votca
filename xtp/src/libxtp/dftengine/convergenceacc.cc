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

#include <algorithm>

// Local VOTCA includes
#include "votca/xtp/convergenceacc.h"
#include "votca/xtp/screening_kernels.h"

namespace votca {
namespace xtp {

/**
 * Convergence acceleration for the SCF cycle.
 *
 * The central residual is the AO-metric commutator R = F P S - S P F, which
 * vanishes at self-consistency in a non-orthogonal basis. This file provides
 * the orthogonalized residual, DIIS or ADIIS extrapolation, level shifting,
 * and occupation-model-specific density construction.
 */

// Build the symmetric orthogonalization matrix X = S^{-1/2}. All Fock-like
// matrices are diagonalized in the orthogonal AO basis X^T F X.
void ConvergenceAcc::setOverlap(AOOverlap& S, double etol) {
  S_ = &S;
  // Own eigendecomposition rather than AOOverlap::Pseudo_InvSqrt: the
  // removed directions are needed as well, see removed_projector_.
  const SymmetricEigenSystem es = SymmetricEigen(S.Matrix());
  const Eigen::VectorXd& s_eig = es.values;
  Eigen::VectorXd inv_sqrt = Eigen::VectorXd::Zero(s_eig.size());
  Index removed = 0;
  for (Index i = 0; i < s_eig.size(); ++i) {
    if (s_eig(i) < etol) {
      ++removed;
    } else {
      inv_sqrt(i) = 1.0 / std::sqrt(s_eig(i));
    }
  }
  Sminusahalf = es.vectors * inv_sqrt.asDiagonal() * es.vectors.transpose();
  removed_projector_.resize(0, 0);
  if (removed > 0) {
    // In SolveFockmatrix, H_ortho = X^T H X vanishes on these directions,
    // so they would come out as eigenvalue-0 "orbitals" in the middle of
    // the spectrum. Shifting them far up keeps them out of the occupied
    // and low virtual space.
    const Eigen::MatrixXd U = es.vectors.leftCols(removed);
    removed_projector_ = U * U.transpose();
  }
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Smallest value of AOOverlap matrix is " << s_eig(0)
      << std::flush;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Removed " << removed
      << " basisfunction from inverse overlap matrix (threshold " << etol << ")"
      << std::flush;
  if (removed > 0) {
    XTP_LOG(Log::error, *log_)
        << TimeStamp() << " The " << removed
        << " removed direction(s) appear as zero orbitals at +" << kRemovedShift
        << " Hartree." << std::flush;
  } else if (s_eig(0) < 1e3 * etol) {
    XTP_LOG(Log::error, *log_)
        << TimeStamp()
        << " WARNING: the overlap matrix is nearly singular; if the SCF is "
           "unstable, raise xtpdft.overlap_tolerance above "
        << s_eig(0) << "." << std::flush;
  }
  return;
}

// Perform one SCF acceleration step.
//
// The commutator residual
//
//   R = X^T (F P S - S P F) X
//
// vanishes at self-consistency in a non-orthogonal AO basis. Its maximum
// element is used as the DIIS error metric, while the history of Fock and
// density matrices is passed to either DIIS or ADIIS to construct the next
// extrapolated Fock matrix. When the extrapolation is deemed unsafe, linear
// density mixing is used instead.
Eigen::MatrixXd ConvergenceAcc::Iterate(const Eigen::MatrixXd& dmat,
                                        Eigen::MatrixXd& H,
                                        tools::EigenSystem& MOs, double totE) {
  totE_.push_back(totE);
  // The MOs of the previous step, for the level shift.
  const Eigen::MatrixXd MOs_old = MOs.eigenvectors();
  const Eigen::VectorXd MOs_old_energies = MOs.eigenvalues();

  const Eigen::MatrixXd& S = S_->Matrix();
  const Eigen::MatrixXd errormatrix =
      Sminusahalf.transpose() * (H * dmat * S - S * dmat * H) * Sminusahalf;
  diiserror_ = errormatrix.cwiseAbs().maxCoeff();

  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " DIIs error " << getDIIsError() << std::flush;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Delta Etot " << getDeltaE() << std::flush;

  // Energy-rise reset. An extrapolated Fock matrix that sends the energy
  // far above anything seen before is not a step to build on: the history
  // that produced it is discarded and the SCF continues, damped, from the
  // lowest-energy density so far.
  ++iterations_since_reset_;
  if (opt_.energy_reset > 0.0 && have_best_ &&
      totE > best_energy_ + opt_.energy_reset &&
      iterations_since_reset_ > kResetCooldown &&
      energy_resets_ < kMaxEnergyResets) {
    ++energy_resets_;
    iterations_since_reset_ = 0;
    XTP_LOG(Log::error, *log_)
        << TimeStamp() << " WARNING: energy rose by " << totE - best_energy_
        << " Ha above the lowest so far (" << best_energy_
        << " Ha); discarding the (A)DIIS history and restarting from that "
           "density (reset "
        << energy_resets_ << ")" << std::flush;
    mathist_.clear();
    dmatHist_.clear();
    errhist_.clear();
    diis_.Clear();
    Eigen::MatrixXd H_restart = best_H_;
    if (opt_.levelshift > 0.0) {
      Levelshift(H_restart, MOs_old);
    }
    MOs = SolveFockmatrix(H_restart);
    usedmixing_ = true;
    return opt_.mixingparameter * best_dmat_ +
           (1.0 - opt_.mixingparameter) * DensityMatrix(MOs);
  }
  if (!have_best_ || totE < best_energy_) {
    have_best_ = true;
    best_energy_ = totE;
    best_dmat_ = dmat;
    best_H_ = H;
  }

  // History, trimmed at one index for all of mathist_, dmatHist_, errhist_
  // and DIIS's own error history: the oldest entry, or with DIIS_maxout the
  // one with the largest error. All four are always the same length, so
  // DIIS::Update trims at the same moment and the same index.
  Index drop = 0;
  if (Index(mathist_.size()) == opt_.histlength) {
    if (opt_.maxout) {
      drop = Index(std::max_element(errhist_.begin(), errhist_.end()) -
                   errhist_.begin());
    }
    mathist_.erase(mathist_.begin() + drop);
    dmatHist_.erase(dmatHist_.begin() + drop);
    errhist_.erase(errhist_.begin() + drop);
  }
  // The UNSHIFTED Fock matrix goes into the history: the level shift is a
  // device for the next diagonalization only, and would otherwise enter
  // the DIIS errors and the ADIIS energy model.
  mathist_.push_back(H);
  dmatHist_.push_back(dmat);
  errhist_.push_back(diiserror_);
  diis_.Update(drop, errormatrix);

  bool diis_error = false;
  Eigen::MatrixXd H_guess = H;
  // Below DIIS_start the density is already close to self-consistency (as
  // after a warm start from the previous QM/MM iteration): no damping, a
  // plain step from the current Fock matrix first and DIIS from two
  // history entries on. Far from it (a cold start), the first steps are
  // damped until the history holds three entries, as before.
  const bool near_convergence = diiserror_ < opt_.diis_start;
  const std::size_t min_history = near_convergence ? 1 : 2;
  if ((diiserror_ < opt_.adiis_start || diiserror_ < opt_.diis_start) &&
      opt_.usediis && mathist_.size() > min_history) {
    Eigen::VectorXd coeffs;
    // ADIIS above DIIS_start, and also below it whenever the energy went
    // up: plain DIIS does not minimize the energy and can run away from a
    // bad step.
    const bool energy_rose = getDeltaE() > kEnergyRiseForADIIS;
    if (diiserror_ > opt_.diis_start || energy_rose) {
      coeffs = adiis_.CalcCoeff(dmatHist_, mathist_);
      diis_error = !adiis_.Info();
      XTP_LOG(Log::warning, *log_)
          << TimeStamp() << " Using ADIIS for next guess" << std::flush;
    } else {
      coeffs = diis_.CalcCoeff();
      diis_error = !diis_.Info();
      XTP_LOG(Log::warning, *log_)
          << TimeStamp() << " Using DIIS for next guess" << std::flush;
    }
    if (diis_error) {
      XTP_LOG(Log::warning, *log_)
          << TimeStamp() << " (A)DIIS failed using mixing instead"
          << std::flush;
    } else {
      H_guess.setZero();
      for (Index i = 0; i < coeffs.size(); i++) {
        if (std::abs(coeffs(i)) < 1e-8) {
          continue;
        }
        H_guess += coeffs(i) * mathist_[i];
      }
    }
  }

  if (opt_.mode != KSmode::fractional && nocclevels_ > 0 &&
      nocclevels_ < MOs_old_energies.size()) {
    const double gap =
        MOs_old_energies(nocclevels_) - MOs_old_energies(nocclevels_ - 1);
    if ((diiserror_ > opt_.levelshiftend && opt_.levelshift > 0.0) ||
        gap < 1e-6) {
      Levelshift(H_guess, MOs_old);
    }
  }

  MOs = SolveFockmatrix(H_guess);

  Eigen::MatrixXd dmatout = DensityMatrix(MOs);
  // mixing_end, not ADIIS_start, decides on damping (as in the UKS path):
  // the two are separate options.
  if (diiserror_ > opt_.mixingend || !opt_.usediis || diis_error ||
      (mathist_.size() <= 2 && !near_convergence)) {
    usedmixing_ = true;
    dmatout =
        opt_.mixingparameter * dmat + (1.0 - opt_.mixingparameter) * dmatout;
    XTP_LOG(Log::warning, *log_)
        << TimeStamp() << " Using Mixing with alpha=" << opt_.mixingparameter
        << std::flush;
  } else {
    usedmixing_ = false;
  }
  return dmatout;
}

bool ConvergenceAcc::HistoryIsAligned() const {
  const std::vector<Eigen::MatrixXd>& errors = diis_.ErrorHistory();
  if (errors.size() != mathist_.size() || dmatHist_.size() != mathist_.size() ||
      errhist_.size() != mathist_.size()) {
    return false;
  }
  const Eigen::MatrixXd& S = S_->Matrix();
  for (std::size_t k = 0; k < mathist_.size(); ++k) {
    const Eigen::MatrixXd e =
        Sminusahalf.transpose() *
        (mathist_[k] * dmatHist_[k] * S - S * dmatHist_[k] * mathist_[k]) *
        Sminusahalf;
    if ((e - errors[k]).cwiseAbs().maxCoeff() > 1e-12 ||
        std::abs(e.cwiseAbs().maxCoeff() - errhist_[k]) > 1e-12) {
      return false;
    }
  }
  return true;
}

void ConvergenceAcc::PrintConfigOptions() const {
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Convergence Options:" << std::flush;
  XTP_LOG(Log::error, *log_)
      << "\t\t Delta E [Ha]: " << opt_.Econverged << std::flush;
  XTP_LOG(Log::error, *log_)
      << "\t\t DIIS max error: " << opt_.error_converged << std::flush;
  if (opt_.usediis) {
    XTP_LOG(Log::error, *log_)
        << "\t\t DIIS histlength: " << opt_.histlength << std::flush;
    XTP_LOG(Log::error, *log_)
        << "\t\t ADIIS start: " << opt_.adiis_start << std::flush;
    XTP_LOG(Log::error, *log_)
        << "\t\t DIIS start: " << opt_.diis_start << std::flush;
    std::string del = "oldest";
    if (opt_.maxout) {
      del = "largest";
    }
    XTP_LOG(Log::error, *log_)
        << "\t\t Deleting " << del << " element from DIIS hist" << std::flush;
  }
  XTP_LOG(Log::error, *log_)
      << "\t\t Levelshift[Ha]: " << opt_.levelshift << std::flush;
  XTP_LOG(Log::error, *log_)
      << "\t\t Levelshift end: " << opt_.levelshiftend << std::flush;
  XTP_LOG(Log::error, *log_)
      << "\t\t Mixing Parameter alpha: " << opt_.mixingparameter << std::flush;
  XTP_LOG(Log::error, *log_)
      << "\t\t Mixing end: " << opt_.mixingend << std::flush;
  XTP_LOG(Log::error, *log_)
      << "\t\t Energy reset [Ha]: " << opt_.energy_reset
      << (opt_.energy_reset > 0.0 ? "" : " (off)") << std::flush;
}

// Solve the generalized Roothaan-Hall problem
//
//   F C = S C eps
//
// by symmetric orthogonalization: diagonalize X^T F X with X = S^{-1/2} and
// back-transform the eigenvectors as C = X C' .
tools::EigenSystem ConvergenceAcc::SolveFockmatrix(
    const Eigen::MatrixXd& H) const {
  // transform to orthogonal for
  Eigen::MatrixXd H_ortho = Sminusahalf.transpose() * H * Sminusahalf;
  if (removed_projector_.size() > 0) {
    H_ortho += kRemovedShift * removed_projector_;
  }
  const SymmetricEigenSystem es = SymmetricEigen(H_ortho);
  if (!es.values.allFinite()) {
    throw std::runtime_error("Matrix Diagonalisation failed");
  }

  tools::EigenSystem result;
  result.eigenvalues() = es.values;
  result.eigenvectors() = Sminusahalf * es.vectors;
  return result;
}

// Add a virtual-space level shift
//
//   F <- F + S C_virt Delta C_virt^T S,
//
// implemented by placing the scalar shift on the virtual diagonal in the MO
// basis built from the previous iteration. This leaves occupied orbitals
// unchanged while opening the HOMO-LUMO gap during difficult SCF phases.
void ConvergenceAcc::Levelshift(Eigen::MatrixXd& H,
                                const Eigen::MatrixXd& MOs_old) const {
  if (opt_.levelshift < 1e-9) {
    return;
  }
  Eigen::VectorXd virt = Eigen::VectorXd::Zero(H.rows());
  for (Index i = nocclevels_; i < H.rows(); i++) {
    virt(i) = opt_.levelshift;
  }

  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Using levelshift:" << opt_.levelshift << " Hartree"
      << std::flush;
  Eigen::MatrixXd vir = S_->Matrix() * MOs_old * virt.asDiagonal() *
                        MOs_old.transpose() * S_->Matrix();
  H += vir;
  return;
}

/*
Eigen::MatrixXd ConvergenceAcc::DensityMatrix(
    const tools::EigenSystem& MOs) const {
  Eigen::MatrixXd result;
  if (opt_.mode == KSmode::closed) {
    result = DensityMatrixGroundState(MOs.eigenvectors());
  } else if (opt_.mode == KSmode::open) {
    result = DensityMatrixGroundState_unres(MOs.eigenvectors());
  } else if (opt_.mode == KSmode::fractional) {
    result = DensityMatrixGroundState_frac(MOs);
  }
  return result;
} */

// Closed-shell AO density matrix
//
//   P = 2 C_occ C_occ^T,
//
// where each occupied spatial orbital contributes one alpha and one beta
// electron.
Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState(
    const Eigen::MatrixXd& MOs) const {
  const Eigen::MatrixXd occstates = MOs.leftCols(nocclevels_);
  Eigen::MatrixXd dmatGS = 2.0 * occstates * occstates.transpose();
  return dmatGS;
}

// Spin-resolved unrestricted AO density matrix for one spin channel,
//
//   P^sigma = C_occ^sigma (C_occ^sigma)^T,
//
// with no factor of two because a single spin channel is represented.
Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState_unres(
    const Eigen::MatrixXd& MOs) const {
  if (nocclevels_ == 0) {
    return Eigen::MatrixXd::Zero(MOs.rows(), MOs.rows());
  }
  Eigen::MatrixXd occstates = MOs.leftCols(nocclevels_);
  Eigen::MatrixXd dmatGS = occstates * occstates.transpose();
  return dmatGS;
}

// Fractionally occupied AO density matrix
//
//   P = C f C^T,
//
// where f is a diagonal matrix of orbital occupations assembled from the
// configured electron count.
// Fractional-occupation AO density matrix
//
//   P = C n C^T,
//
// where n is the diagonal matrix of orbital occupations provided in the
// EigenSystem container.
Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState_frac(
    const tools::EigenSystem& MOs) const {
  if (opt_.numberofelectrons == 0) {
    return Eigen::MatrixXd::Zero(MOs.eigenvectors().rows(),
                                 MOs.eigenvectors().rows());
  }

  Eigen::VectorXd occupation = Eigen::VectorXd::Zero(MOs.eigenvalues().size());
  std::vector<std::vector<Index> > degeneracies;
  double buffer = 1e-4;
  degeneracies.push_back(std::vector<Index>{0});
  for (Index i = 1; i < occupation.size(); i++) {
    if (MOs.eigenvalues()(i) <
        MOs.eigenvalues()(degeneracies[degeneracies.size() - 1][0]) + buffer) {
      degeneracies[degeneracies.size() - 1].push_back(i);
    } else {
      degeneracies.push_back(std::vector<Index>{i});
    }
  }
  Index numofelec = opt_.numberofelectrons;
  for (const std::vector<Index>& deglevel : degeneracies) {
    Index numofpossibleelectrons = 2 * Index(deglevel.size());
    if (numofpossibleelectrons <= numofelec) {
      for (Index i : deglevel) {
        occupation(i) = 2;
      }
      numofelec -= numofpossibleelectrons;
    } else {
      double occ = double(numofelec) / double(deglevel.size());
      for (Index i : deglevel) {
        occupation(i) = occ;
      }
      break;
    }
  }
  Eigen::MatrixXd dmatGS = MOs.eigenvectors() * occupation.asDiagonal() *
                           MOs.eigenvectors().transpose();
  return dmatGS;
}

/*******************************************************
 * EXTENSION FOR SPIN-KS-DFT
 *******************************************************/

ConvergenceAcc::SpinDensity
    ConvergenceAcc::DensityMatrixGroundState_restricted_open(
        const Eigen::MatrixXd& MOs) const {

  const Index n_docc =
      std::min(opt_.number_alpha_electrons, opt_.number_beta_electrons);
  const Index n_socc_alpha = opt_.number_alpha_electrons - n_docc;

  SpinDensity result;
  result.alpha = Eigen::MatrixXd::Zero(MOs.rows(), MOs.rows());
  result.beta = Eigen::MatrixXd::Zero(MOs.rows(), MOs.rows());

  if (n_docc > 0) {
    const Eigen::MatrixXd docc = MOs.leftCols(n_docc);
    const Eigen::MatrixXd d_docc = docc * docc.transpose();
    result.alpha += d_docc;
    result.beta += d_docc;
  }

  if (n_socc_alpha > 0) {
    const Eigen::MatrixXd socc = MOs.middleCols(n_docc, n_socc_alpha);
    result.alpha += socc * socc.transpose();
  }

  return result;
}

// Construct spin-resolved densities according to the configured occupation
// model. For restricted open-shell cases the same spatial orbitals are split
// into doubly and singly occupied subsets using the stored alpha/beta counts.
ConvergenceAcc::SpinDensity ConvergenceAcc::DensityMatrixSpinResolved(
    const tools::EigenSystem& MOs) const {

  if (opt_.mode == KSmode::restricted_open) {
    return DensityMatrixGroundState_restricted_open(MOs.eigenvectors());
  } else if (opt_.mode == KSmode::closed) {
    Eigen::MatrixXd d = DensityMatrixGroundState(MOs.eigenvectors());
    return {0.5 * d, 0.5 * d};
  } else if (opt_.mode == KSmode::open) {
    Eigen::MatrixXd d = DensityMatrixGroundState_unres(MOs.eigenvectors());
    return {d, Eigen::MatrixXd::Zero(d.rows(), d.cols())};
  } else {
    Eigen::MatrixXd d = DensityMatrixGroundState_frac(MOs);
    return {0.5 * d, 0.5 * d};
  }
}

Eigen::MatrixXd ConvergenceAcc::DensityMatrix(
    const tools::EigenSystem& MOs) const {
  SpinDensity spin_dmat = DensityMatrixSpinResolved(MOs);
  return spin_dmat.total();
}

}  // namespace xtp
}  // namespace votca
