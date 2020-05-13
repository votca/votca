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
#include "votca/xtp/adiis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/diis.h"
#include "votca/xtp/logger.h"

namespace votca {
namespace xtp {

class ConvergenceAcc {
 public:
  enum KSmode { closed, open, fractional };

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
  };

  void Configure(const ConvergenceAcc::options& opt) {
    _opt = opt;
    if (_opt.mode == KSmode::closed) {
      _nocclevels = _opt.numberofelectrons / 2;
    } else if (_opt.mode == KSmode::open) {
      _nocclevels = _opt.numberofelectrons;
    } else if (_opt.mode == KSmode::fractional) {
      _nocclevels = 0;
    }
    _diis.setHistLength(_opt.histlength);
  }
  void setLogger(Logger* log) { _log = log; }

  void PrintConfigOptions() const;

  bool isConverged() const {
    if (_totE.size() < 2) {
      return false;
    } else {
      return std::abs(getDeltaE()) < _opt.Econverged &&
             getDIIsError() < _opt.error_converged;
    }
  }

  double getDeltaE() const {
    if (_totE.size() < 2) {
      return 0;
    } else {
      return _totE.back() - _totE[_totE.size() - 2];
    }
  }
  void setOverlap(AOOverlap& S, double etol);

  double getDIIsError() const { return _diiserror; }

  bool getUseMixing() const { return _usedmixing; }

  Eigen::MatrixXd Iterate(const Eigen::MatrixXd& dmat, Eigen::MatrixXd& H,
                          tools::EigenSystem& MOs, double totE);
  tools::EigenSystem SolveFockmatrix(const Eigen::MatrixXd& H) const;
  void Levelshift(Eigen::MatrixXd& H, const Eigen::MatrixXd& MOs_old) const;

  Eigen::MatrixXd DensityMatrix(const tools::EigenSystem& MOs) const;

 private:
  options _opt;

  Eigen::MatrixXd DensityMatrixGroundState(const Eigen::MatrixXd& MOs) const;
  Eigen::MatrixXd DensityMatrixGroundState_unres(
      const Eigen::MatrixXd& MOs) const;
  Eigen::MatrixXd DensityMatrixGroundState_frac(
      const tools::EigenSystem& MOs) const;

  bool _usedmixing = true;
  double _diiserror = std::numeric_limits<double>::max();
  Logger* _log;
  const AOOverlap* _S;

  Eigen::MatrixXd Sminusahalf;
  std::vector<Eigen::MatrixXd> _mathist;
  std::vector<Eigen::MatrixXd> _dmatHist;
  std::vector<double> _totE;

  Index _nocclevels;
  Index _maxerrorindex = 0;
  double _maxerror = 0.0;
  ADIIS _adiis;
  DIIS _diis;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CONVERGENCEACC_H
