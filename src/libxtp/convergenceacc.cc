/*
 *            Copyright 2009-2018 The VOTCA Development Team
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
#include "votca/xtp/convergenceacc.h"
#include "votca/xtp/aomatrix.h"

namespace votca {
namespace xtp {

void ConvergenceAcc::setOverlap(AOOverlap& S, double etol) {
  _S = &S;
  Sminusahalf = S.Pseudo_InvSqrt(etol);
  CTP_LOG(ctp::logDEBUG, *_log)
      << ctp::TimeStamp() << " Smallest value of AOOverlap matrix is "
      << _S->SmallestEigenValue() << std::flush;
  CTP_LOG(ctp::logDEBUG, *_log)
      << ctp::TimeStamp() << " Removed " << _S->Removedfunctions()
      << " basisfunction from inverse overlap matrix" << std::flush;
  Sonehalf = S.Sqrt();
  return;
}

Eigen::MatrixXd ConvergenceAcc::Iterate(const Eigen::MatrixXd& dmat,
                                        Eigen::MatrixXd& H,
                                        Eigen::VectorXd& MOenergies,
                                        Eigen::MatrixXd& MOs, double totE) {
  Eigen::MatrixXd H_guess = Eigen::MatrixXd::Zero(H.rows(), H.cols());

  if (int(_mathist.size()) == _opt.histlength) {
    _totE.erase(_totE.begin() + _maxerrorindex);
    _mathist.erase(_mathist.begin() + _maxerrorindex);
    _dmatHist.erase(_dmatHist.begin() + _maxerrorindex);
  }

  _totE.push_back(totE);
  if (_opt.mode != KSmode::fractional) {
    double gap = MOenergies(_nocclevels) - MOenergies(_nocclevels - 1);
    if ((_diiserror > _opt.levelshiftend && _opt.levelshift > 0.0) ||
        gap < 1e-6) {
      Levelshift(H);
    }
  }
  const Eigen::MatrixXd& S = _S->Matrix();
  Eigen::MatrixXd errormatrix =
      Sminusahalf.transpose() * (H * dmat * S - S * dmat * H) * Sminusahalf;
  _diiserror = errormatrix.cwiseAbs().maxCoeff();

  _mathist.push_back(H);
  _dmatHist.push_back(dmat);

  if (_opt.maxout) {
    if (_diiserror > _maxerror) {
      _maxerror = _diiserror;
      _maxerrorindex = _mathist.size() - 1;
    }
  }

  _diis.Update(_maxerrorindex, errormatrix);
  bool diis_error = false;
  CTP_LOG(ctp::logDEBUG, *_log)
      << ctp::TimeStamp() << " DIIs error " << getDIIsError() << std::flush;

  CTP_LOG(ctp::logDEBUG, *_log)
      << ctp::TimeStamp() << " Delta Etot " << getDeltaE() << std::flush;

  if ((_diiserror < _opt.adiis_start || _diiserror < _opt.diis_start) &&
      _opt.usediis && _mathist.size() > 2) {
    Eigen::VectorXd coeffs;
    // use ADIIs if energy has risen a lot in current iteration

    if (_diiserror > _opt.diis_start ||
        _totE.back() > 0.9 * _totE[_totE.size() - 2]) {
      coeffs = _adiis.CalcCoeff(_dmatHist, _mathist);
      diis_error = !_adiis.Info();
      CTP_LOG(ctp::logDEBUG, *_log)
          << ctp::TimeStamp() << " Using ADIIS for next guess" << std::flush;

    } else {
      coeffs = _diis.CalcCoeff();
      diis_error = !_diis.Info();
      CTP_LOG(ctp::logDEBUG, *_log)
          << ctp::TimeStamp() << " Using DIIS for next guess" << std::flush;
    }
    if (diis_error) {
      CTP_LOG(ctp::logDEBUG, *_log)
          << ctp::TimeStamp() << " (A)DIIS failed using mixing instead"
          << std::flush;
      H_guess = H;
    } else {
      for (int i = 0; i < coeffs.size(); i++) {
        if (std::abs(coeffs(i)) < 1e-8) {
          continue;
        }
        H_guess += coeffs(i) * _mathist[i];
      }
    }

  } else {
    H_guess = H;
  }

  SolveFockmatrix(MOenergies, MOs, H_guess);
  Eigen::MatrixXd dmatout = DensityMatrix(MOs, MOenergies);

  if (_diiserror > _opt.adiis_start || !_opt.usediis || diis_error ||
      _mathist.size() <= 2) {
    _usedmixing = true;
    dmatout =
        _opt.mixingparameter * dmat + (1.0 - _opt.mixingparameter) * dmatout;
    CTP_LOG(ctp::logDEBUG, *_log)
        << ctp::TimeStamp()
        << " Using Mixing with alpha=" << _opt.mixingparameter << std::flush;
  } else {
    _usedmixing = false;
  }
  return dmatout;
}

void ConvergenceAcc::PrintConfigOptions() const {
  CTP_LOG(ctp::logDEBUG, *_log)
      << ctp::TimeStamp() << " Convergence Options:" << std::flush;
  CTP_LOG(ctp::logDEBUG, *_log)
      << "\t\t Delta E [Ha]: " << _opt.Econverged << std::flush;
  CTP_LOG(ctp::logDEBUG, *_log)
      << "\t\t DIIS max error: " << _opt.error_converged << std::flush;
  if (_opt.usediis) {
    CTP_LOG(ctp::logDEBUG, *_log)
        << "\t\t DIIS histlength: " << _opt.histlength << std::flush;
    CTP_LOG(ctp::logDEBUG, *_log)
        << "\t\t ADIIS start: " << _opt.adiis_start << std::flush;
    CTP_LOG(ctp::logDEBUG, *_log)
        << "\t\t DIIS start: " << _opt.diis_start << std::flush;
    std::string del = "oldest";
    if (_opt.maxout) {
      del = "largest";
    }
    CTP_LOG(ctp::logDEBUG, *_log)
        << "\t\t Deleting " << del << " element from DIIS hist" << std::flush;
  }
  CTP_LOG(ctp::logDEBUG, *_log)
      << "\t\t Levelshift[Ha]: " << _opt.levelshift << std::flush;
  CTP_LOG(ctp::logDEBUG, *_log)
      << "\t\t Levelshift end: " << _opt.levelshiftend << std::flush;
  CTP_LOG(ctp::logDEBUG, *_log)
      << "\t\t Mixing Parameter alpha: " << _opt.mixingparameter << std::flush;
}

void ConvergenceAcc::SolveFockmatrix(Eigen::VectorXd& MOenergies,
                                     Eigen::MatrixXd& MOs,
                                     const Eigen::MatrixXd& H) {
  // transform to orthogonal for
  Eigen::MatrixXd H_ortho = Sminusahalf.transpose() * H * Sminusahalf;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_ortho);

  if (es.info() != Eigen::ComputationInfo::Success) {
    throw std::runtime_error("Matrix Diagonalisation failed. DiagInfo" +
                             std::to_string(es.info()));
  }

  MOsinv = es.eigenvectors().transpose() * Sonehalf;
  MOenergies = es.eigenvalues();

  MOs = Sminusahalf * es.eigenvectors();
  return;
}

void ConvergenceAcc::Levelshift(Eigen::MatrixXd& H) const {
  if (_opt.levelshift < 1e-9) {
    return;
  }
  if (MOsinv.rows() < 1) {
    throw std::runtime_error(
        "ConvergenceAcc::Levelshift: Call SolveFockmatrix before Levelshift, "
        "MOsinv not initialized");
  }
  Eigen::MatrixXd virt = Eigen::MatrixXd::Zero(H.rows(), H.cols());
  for (int i = _nocclevels; i < H.rows(); i++) {
    virt(i, i) = _opt.levelshift;
  }

  CTP_LOG(ctp::logDEBUG, *_log)
      << ctp::TimeStamp() << " Using levelshift:" << _opt.levelshift
      << " Hartree" << std::flush;
  Eigen::MatrixXd vir = MOsinv.transpose() * virt * MOsinv;
  H += vir;
  return;
}

Eigen::MatrixXd ConvergenceAcc::DensityMatrix(
    const Eigen::MatrixXd& MOs, const Eigen::VectorXd& MOEnergies) const {
  Eigen::MatrixXd result;
  if (_opt.mode == KSmode::closed) {
    result = DensityMatrixGroundState(MOs);
  } else if (_opt.mode == KSmode::open) {
    result = DensityMatrixGroundState_unres(MOs);
  } else if (_opt.mode == KSmode::fractional) {
    result = DensityMatrixGroundState_frac(MOs, MOEnergies);
  }
  return result;
}

Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState(
    const Eigen::MatrixXd& MOs) const {
  const Eigen::MatrixXd occstates = MOs.block(0, 0, MOs.rows(), _nocclevels);
  Eigen::MatrixXd dmatGS = 2.0 * occstates * occstates.transpose();
  return dmatGS;
}

Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState_unres(
    const Eigen::MatrixXd& MOs) const {
  if (_nocclevels == 0) {
    return Eigen::MatrixXd::Zero(MOs.cols(), MOs.rows());
  }
  Eigen::MatrixXd occstates = MOs.block(0, 0, MOs.rows(), _nocclevels);
  Eigen::MatrixXd dmatGS = occstates * occstates.transpose();
  return dmatGS;
}

Eigen::MatrixXd ConvergenceAcc::DensityMatrixGroundState_frac(
    const Eigen::MatrixXd& MOs, const Eigen::VectorXd& MOEnergies) const {
  if (_opt.numberofelectrons == 0) {
    return Eigen::MatrixXd::Zero(MOs.rows(), MOs.cols());
  }
  int numofelec = _opt.numberofelectrons;
  Eigen::VectorXd occupation = Eigen::VectorXd::Zero(MOEnergies.size());

  std::vector<std::vector<int> > degeneracies;
  double buffer = 1e-4;
  degeneracies.push_back(std::vector<int>{0});
  for (int i = 1; i < occupation.size(); i++) {
    if (MOEnergies(i) <
        MOEnergies(degeneracies[degeneracies.size() - 1][0]) + buffer) {
      degeneracies[degeneracies.size() - 1].push_back(i);
    } else {
      degeneracies.push_back(std::vector<int>{i});
    }
  }
  for (const std::vector<int>& deglevel : degeneracies) {
    int numofpossibleelectrons = 2 * deglevel.size();
    if (numofpossibleelectrons <= numofelec) {
      for (const int& i : deglevel) {
        occupation(i) = 2;
      }
      numofelec -= numofpossibleelectrons;
    } else if (numofpossibleelectrons > numofelec) {
      double occ = double(numofelec) / double(deglevel.size());
      for (const int& i : deglevel) {
        occupation(i) = occ;
      }
      break;
    }
  }
  Eigen::MatrixXd dmatGS = MOs * occupation.asDiagonal() * MOs.transpose();
  return dmatGS;
}

}  // namespace xtp
}  // namespace votca
