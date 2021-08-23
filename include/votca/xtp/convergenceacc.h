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
    opt_ = opt;
    if (opt_.mode == KSmode::closed) {
      nocclevels_ = opt_.numberofelectrons / 2;
    } else if (opt_.mode == KSmode::open) {
      nocclevels_ = opt_.numberofelectrons;
    } else if (opt_.mode == KSmode::fractional) {
      nocclevels_ = 0;
    }
    diis_.setHistLength(opt_.histlength);
  }
  void setLogger(Logger* log) { log_ = log; }

  void PrintConfigOptions() const;

  bool isConverged() const {
    if (totE_.size() < 2) {
      return false;
    } else {
      return std::abs(getDeltaE()) < opt_.Econverged &&
             getDIIsError() < opt_.error_converged;
    }
  }

  double getDeltaE() const {
    if (totE_.size() < 2) {
      return 0;
    } else {
      return totE_.back() - totE_[totE_.size() - 2];
    }
  }
  void setOverlap(AOOverlap& S, double etol);

  double getDIIsError() const { return diiserror_; }

  bool getUseMixing() const { return usedmixing_; }

  Eigen::MatrixXd Iterate(const Eigen::MatrixXd& dmat, Eigen::MatrixXd& H,
                          tools::EigenSystem& MOs, double totE);
  tools::EigenSystem SolveFockmatrix(const Eigen::MatrixXd& H) const;
  void Levelshift(Eigen::MatrixXd& H, const Eigen::MatrixXd& MOs_old) const;

  Eigen::MatrixXd DensityMatrix(const tools::EigenSystem& MOs) const;

 private:
  options opt_;

  Eigen::MatrixXd DensityMatrixGroundState(const Eigen::MatrixXd& MOs) const;
  Eigen::MatrixXd DensityMatrixGroundState_unres(
      const Eigen::MatrixXd& MOs) const;
  Eigen::MatrixXd DensityMatrixGroundState_frac(
      const tools::EigenSystem& MOs) const;

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
