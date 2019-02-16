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

#ifndef _VOTCA_XTP_CONVERGENCEACC__H
#define _VOTCA_XTP_CONVERGENCEACC__H

#include <memory>
#include <votca/ctp/logger.h>
#include <votca/xtp/adiis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/basisset.h>
#include <votca/xtp/diis.h>
namespace votca {
namespace xtp {

class ConvergenceAcc {
 public:
  enum KSmode { closed, open, fractional };

  struct options {
    KSmode mode = KSmode::closed;
    bool usediis = true;
    bool noisy = false;
    int histlength = 10;
    bool maxout = false;
    double adiis_start = 2;
    double diis_start = 0.01;
    double levelshift = 0.25;
    double levelshiftend = 0.8;
    int numberofelectrons;
    double mixingparameter = 0.7;
    double Econverged = 1e-7;
    double error_converged = 1e-7;
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
  void setLogger(ctp::Logger* log) { _log = log; }

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
                          Eigen::VectorXd& MOenergies, Eigen::MatrixXd& MOs,
                          double totE);
  void SolveFockmatrix(Eigen::VectorXd& MOenergies, Eigen::MatrixXd& MOs,
                       const Eigen::MatrixXd& H);
  void Levelshift(Eigen::MatrixXd& H) const;

  Eigen::MatrixXd DensityMatrix(const Eigen::MatrixXd& MOs,
                                const Eigen::VectorXd& MOEnergies) const;

 private:
  options _opt;

  Eigen::MatrixXd DensityMatrixGroundState(const Eigen::MatrixXd& MOs) const;
  Eigen::MatrixXd DensityMatrixGroundState_unres(
      const Eigen::MatrixXd& MOs) const;
  Eigen::MatrixXd DensityMatrixGroundState_frac(
      const Eigen::MatrixXd& MOs, const Eigen::VectorXd& MOEnergies) const;

  bool _usedmixing = true;
  double _diiserror = std::numeric_limits<double>::max();
  ctp::Logger* _log;
  const AOOverlap* _S;

  Eigen::MatrixXd Sminusahalf;
  Eigen::MatrixXd Sonehalf;
  Eigen::MatrixXd MOsinv;
  std::vector<Eigen::MatrixXd> _mathist;
  std::vector<Eigen::MatrixXd> _dmatHist;
  std::vector<double> _totE;

  int _nocclevels;
  int _maxerrorindex = 0;
  double _maxerror = 0.0;
  ADIIS _adiis;
  DIIS _diis;
};

}  // namespace xtp
}  // namespace votca

#endif
