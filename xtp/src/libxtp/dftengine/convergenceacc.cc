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

// Local VOTCA includes
#include "votca/xtp/convergenceacc.h"

namespace votca {
namespace xtp {

void ConvergenceAcc::setOverlap(AOOverlap& S, double etol) {
  S_ = &S;
  Sminusahalf = S.Pseudo_InvSqrt(etol);
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Smallest value of AOOverlap matrix is "
      << S_->SmallestEigenValue() << std::flush;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Removed " << S_->Removedfunctions()
      << " basisfunction from inverse overlap matrix" << std::flush;
  return;
}

Eigen::MatrixXd ConvergenceAcc::Iterate(const Eigen::MatrixXd& dmat,
                                        Eigen::MatrixXd& H,
                                        tools::EigenSystem& MOs, double totE) {
  Eigen::MatrixXd H_guess = Eigen::MatrixXd::Zero(H.rows(), H.cols());

  if (int(mathist_.size()) == opt_.histlength) {
    totE_.erase(totE_.begin() + maxerrorindex_);
    mathist_.erase(mathist_.begin() + maxerrorindex_);
    dmatHist_.erase(dmatHist_.begin() + maxerrorindex_);
  }

  totE_.push_back(totE);
  if (opt_.mode != KSmode::fractional) {
    double gap =
        MOs.eigenvalues()(nocclevels_) - MOs.eigenvalues()(nocclevels_ - 1);
    if ((diiserror_ > opt_.levelshiftend && opt_.levelshift > 0.0) ||
        gap < 1e-6) {
      Levelshift(H, MOs.eigenvectors());
    }
  }
  const Eigen::MatrixXd& S = S_->Matrix();
  Eigen::MatrixXd errormatrix =
      Sminusahalf.transpose() * (H * dmat * S - S * dmat * H) * Sminusahalf;
  diiserror_ = errormatrix.cwiseAbs().maxCoeff();

  mathist_.push_back(H);
  dmatHist_.push_back(dmat);

  if (opt_.maxout) {
    if (diiserror_ > maxerror_) {
      maxerror_ = diiserror_;
      maxerrorindex_ = mathist_.size() - 1;
    }
  }

  diis_.Update(maxerrorindex_, errormatrix);
  bool diis_error = false;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " DIIs error " << getDIIsError() << std::flush;

  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Delta Etot " << getDeltaE() << std::flush;

  if ((diiserror_ < opt_.adiis_start || diiserror_ < opt_.diis_start) &&
      opt_.usediis && mathist_.size() > 2) {
    Eigen::VectorXd coeffs;
    // use ADIIs if energy has risen a lot in current iteration

    if (diiserror_ > opt_.diis_start ||
        totE_.back() > 0.9 * totE_[totE_.size() - 2]) {
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
      H_guess = H;
    } else {
      for (Index i = 0; i < coeffs.size(); i++) {
        if (std::abs(coeffs(i)) < 1e-8) {
          continue;
        }
        H_guess += coeffs(i) * mathist_[i];
      }
    }

  } else {
    H_guess = H;
  }

  MOs = SolveFockmatrix(H_guess);

  // if user selects active{
  // list_activeatoms = pass the list;
  // for (i=0, i<len(list_activeatoms), i++){
  // Eigen::MatrixXd activemos.col(i) = MOs.col(list_activeatoms[i]);
  // }
  // Eigen::MatrixXd dmatout = activemos * activemos.transpose();
  // }
  // else {
  Eigen::MatrixXd dmatout = DensityMatrix(MOs);
  //}
  if (diiserror_ > opt_.adiis_start || !opt_.usediis || diis_error ||
      mathist_.size() <= 2) {
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
}

tools::EigenSystem ConvergenceAcc::SolveFockmatrix(
    const Eigen::MatrixXd& H) const {
  // transform to orthogonal for
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
}

Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState(
    const Eigen::MatrixXd& MOs) const {
  const Eigen::MatrixXd occstates = MOs.leftCols(nocclevels_);
  Eigen::MatrixXd dmatGS = 2.0 * occstates * occstates.transpose();
  return dmatGS;
}

Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState_unres(
    const Eigen::MatrixXd& MOs) const {
  if (nocclevels_ == 0) {
    return Eigen::MatrixXd::Zero(MOs.rows(), MOs.rows());
  }
  Eigen::MatrixXd occstates = MOs.leftCols(nocclevels_);
  Eigen::MatrixXd dmatGS = occstates * occstates.transpose();
  return dmatGS;
}

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

}  // namespace xtp
}  // namespace votca
