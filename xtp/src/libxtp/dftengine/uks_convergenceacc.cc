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

// Local VOTCA includes
#include "votca/xtp/uks_convergenceacc.h"

namespace votca {
namespace xtp {

void UKSConvergenceAcc::Configure(const options& opt_alpha,
                                  const options& opt_beta) {
  opt_alpha_ = opt_alpha;
  opt_beta_ = opt_beta;

  nocclevels_alpha_ = opt_alpha_.numberofelectrons;
  nocclevels_beta_ = opt_beta_.numberofelectrons;

  // one shared DIIS/ADIIS history length
  diis_.setHistLength(opt_alpha_.histlength);
  //adiis_.setHistLength(opt_alpha_.histlength);
}

void UKSConvergenceAcc::setLogger(Logger* log) { log_ = log; }

void UKSConvergenceAcc::setOverlap(AOOverlap& S, double etol) {
  S_ = &S;
  Sminusahalf = S.Pseudo_InvSqrt(etol);
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Smallest value of AOOverlap matrix is "
      << S_->SmallestEigenValue() << std::flush;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Removed " << S_->Removedfunctions()
      << " basisfunction from inverse overlap matrix" << std::flush;
}

tools::EigenSystem UKSConvergenceAcc::SolveFockmatrix(
    const Eigen::MatrixXd& H) const {
  Eigen::MatrixXd H_ortho = Sminusahalf.transpose() * H * Sminusahalf;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_ortho);

  if (es.info() != Eigen::ComputationInfo::Success) {
    throw std::runtime_error("Matrix Diagonalisation failed. DiagInfo" +
                             std::to_string(es.info()));
  }

  tools::EigenSystem result;
  result.eigenvalues() = es.eigenvalues();
  result.eigenvectors() = Sminusahalf * es.eigenvectors();
  return result;
}

Eigen::MatrixXd UKSConvergenceAcc::DensityMatrixGroundState_unres(
    const Eigen::MatrixXd& MOs, Index nocclevels) const {
  if (nocclevels == 0) {
    return Eigen::MatrixXd::Zero(MOs.rows(), MOs.rows());
  }
  Eigen::MatrixXd occstates = MOs.leftCols(nocclevels);
  return occstates * occstates.transpose();
}

UKSConvergenceAcc::SpinDensity UKSConvergenceAcc::DensityMatrix(
    const tools::EigenSystem& MOs_alpha,
    const tools::EigenSystem& MOs_beta) const {
  SpinDensity result;
  result.alpha =
      DensityMatrixGroundState_unres(MOs_alpha.eigenvectors(), nocclevels_alpha_);
  result.beta =
      DensityMatrixGroundState_unres(MOs_beta.eigenvectors(), nocclevels_beta_);
  return result;
}

void UKSConvergenceAcc::Levelshift(Eigen::MatrixXd& H,
                                   const Eigen::MatrixXd& MOs_old,
                                   const options& opt,
                                   Index nocclevels) const {
  if (opt.levelshift < 1e-9) {
    return;
  }
  Eigen::VectorXd virt = Eigen::VectorXd::Zero(H.rows());
  for (Index i = nocclevels; i < H.rows(); ++i) {
    virt(i) = opt.levelshift;
  }

  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Using levelshift:" << opt.levelshift << " Hartree"
      << std::flush;

  Eigen::MatrixXd vir = S_->Matrix() * MOs_old * virt.asDiagonal() *
                        MOs_old.transpose() * S_->Matrix();
  H += vir;
}

Eigen::MatrixXd UKSConvergenceAcc::BuildErrorMatrix(
    const Eigen::MatrixXd& dmat, const Eigen::MatrixXd& H) const {
  const Eigen::MatrixXd& S = S_->Matrix();
  return Sminusahalf.transpose() * (H * dmat * S - S * dmat * H) * Sminusahalf;
}

double UKSConvergenceAcc::CombinedError(const Eigen::MatrixXd& err_alpha,
                                        const Eigen::MatrixXd& err_beta) const {
  return std::max(err_alpha.cwiseAbs().maxCoeff(),
                  err_beta.cwiseAbs().maxCoeff());
}

UKSConvergenceAcc::SpinDensity UKSConvergenceAcc::Iterate(
    const SpinDensity& dmat, SpinFock& H, tools::EigenSystem& MOs_alpha,
    tools::EigenSystem& MOs_beta, double totE) {

  if (int(mathist_alpha_.size()) == opt_alpha_.histlength) {
    totE_.erase(totE_.begin() + maxerrorindex_);
    mathist_alpha_.erase(mathist_alpha_.begin() + maxerrorindex_);
    mathist_beta_.erase(mathist_beta_.begin() + maxerrorindex_);
    dmatHist_alpha_.erase(dmatHist_alpha_.begin() + maxerrorindex_);
    dmatHist_beta_.erase(dmatHist_beta_.begin() + maxerrorindex_);
  }

  totE_.push_back(totE);

  if (nocclevels_alpha_ > 0 && nocclevels_alpha_ < MOs_alpha.eigenvalues().size()) {
    double gap_alpha =
        MOs_alpha.eigenvalues()(nocclevels_alpha_) -
        MOs_alpha.eigenvalues()(nocclevels_alpha_ - 1);
    if ((diiserror_ > opt_alpha_.levelshiftend && opt_alpha_.levelshift > 0.0) ||
        gap_alpha < 1e-6) {
      Levelshift(H.alpha, MOs_alpha.eigenvectors(), opt_alpha_, nocclevels_alpha_);
    }
  }

  if (nocclevels_beta_ > 0 && nocclevels_beta_ < MOs_beta.eigenvalues().size()) {
    double gap_beta =
        MOs_beta.eigenvalues()(nocclevels_beta_) -
        MOs_beta.eigenvalues()(nocclevels_beta_ - 1);
    if ((diiserror_ > opt_beta_.levelshiftend && opt_beta_.levelshift > 0.0) ||
        gap_beta < 1e-6) {
      Levelshift(H.beta, MOs_beta.eigenvectors(), opt_beta_, nocclevels_beta_);
    }
  }

  Eigen::MatrixXd err_alpha = BuildErrorMatrix(dmat.alpha, H.alpha);
  Eigen::MatrixXd err_beta = BuildErrorMatrix(dmat.beta, H.beta);

  diiserror_ = CombinedError(err_alpha, err_beta);

  mathist_alpha_.push_back(H.alpha);
  mathist_beta_.push_back(H.beta);
  dmatHist_alpha_.push_back(dmat.alpha);
  dmatHist_beta_.push_back(dmat.beta);

  if (opt_alpha_.maxout) {
    if (diiserror_ > maxerror_) {
      maxerror_ = diiserror_;
      maxerrorindex_ = mathist_alpha_.size() - 1;
    }
  } else {
    maxerrorindex_ = 0;
  }

  // crucial: one shared error matrix = alpha + beta contribution
  diis_.Update(maxerrorindex_, err_alpha, err_beta);

  bool diis_error = false;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " DIIs error " << diiserror_ << std::flush;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Delta Etot " << getDeltaE() << std::flush;

  Eigen::MatrixXd H_guess_alpha = H.alpha;
  Eigen::MatrixXd H_guess_beta = H.beta;

  if ((diiserror_ < opt_alpha_.adiis_start || diiserror_ < opt_alpha_.diis_start) &&
      opt_alpha_.usediis && mathist_alpha_.size() > 2) {

    Eigen::VectorXd coeffs;

    if (diiserror_ > opt_alpha_.diis_start ||
        totE_.back() > 0.9 * totE_[totE_.size() - 2]) {
      coeffs = adiis_.CalcCoeff(dmatHist_alpha_, dmatHist_beta_,
                                mathist_alpha_, mathist_beta_);
      diis_error = !adiis_.Info() || coeffs.size() == 0;
      XTP_LOG(Log::warning, *log_)
          << TimeStamp() << " Using ADIIS for next UKS guess" << std::flush;
    } else {
coeffs = diis_.CalcCoeff();
diis_error = !diis_.Info() || coeffs.size() == 0;
      XTP_LOG(Log::warning, *log_)
          << TimeStamp() << " Using DIIS for next UKS guess" << std::flush;
    }

    if (diis_error) {
      XTP_LOG(Log::warning, *log_)
          << TimeStamp() << " (A)DIIS failed using mixing instead"
          << std::flush;
      H_guess_alpha = H.alpha;
      H_guess_beta = H.beta;
    } else {
      H_guess_alpha.setZero();
      H_guess_beta.setZero();
      for (Index i = 0; i < coeffs.size(); ++i) {
        if (std::abs(coeffs(i)) < 1e-8) {
          continue;
        }
        H_guess_alpha += coeffs(i) * mathist_alpha_[i];
        H_guess_beta += coeffs(i) * mathist_beta_[i];
      }
    }
  }

  MOs_alpha = SolveFockmatrix(H_guess_alpha);
  MOs_beta = SolveFockmatrix(H_guess_beta);

  SpinDensity dmatout = DensityMatrix(MOs_alpha, MOs_beta);

  if (diiserror_ > opt_alpha_.adiis_start || !opt_alpha_.usediis || diis_error ||
      mathist_alpha_.size() <= 2) {
    usedmixing_ = true;
    dmatout.alpha = opt_alpha_.mixingparameter * dmat.alpha +
                    (1.0 - opt_alpha_.mixingparameter) * dmatout.alpha;
    dmatout.beta = opt_beta_.mixingparameter * dmat.beta +
                   (1.0 - opt_beta_.mixingparameter) * dmatout.beta;
    XTP_LOG(Log::warning, *log_)
        << TimeStamp() << " Using coupled UKS mixing with alpha="
        << opt_alpha_.mixingparameter << std::flush;
  } else {
    usedmixing_ = false;
  }

  return dmatout;
}

double UKSConvergenceAcc::getDeltaE() const {
  if (totE_.size() < 2) {
    return 100.0;
  }
  return std::abs(totE_.back() - totE_[totE_.size() - 2]);
}

bool UKSConvergenceAcc::isConverged() const {
  return (getDeltaE() < opt_alpha_.Econverged &&
          diiserror_ < opt_alpha_.error_converged);
}

}  // namespace xtp
}  // namespace votca