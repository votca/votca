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

#include <vector>
#include <Eigen/Dense>

#include "votca/xtp/aomatrix.h"
#include "votca/xtp/adiis.h"
#include "votca/xtp/convergenceacc.h"
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

  void Configure(const options& opt_alpha, const options& opt_beta);
  void setLogger(Logger* log);
  void setOverlap(AOOverlap& S, double etol);

  SpinDensity DensityMatrix(const tools::EigenSystem& MOs_alpha,
                           const tools::EigenSystem& MOs_beta) const;

  SpinDensity Iterate(const SpinDensity& dmat, SpinFock& H,
                      tools::EigenSystem& MOs_alpha,
                      tools::EigenSystem& MOs_beta, double totE);

  tools::EigenSystem SolveFockmatrix(const Eigen::MatrixXd& H) const;

  bool isConverged() const;
  double getDIIsError() const { return diiserror_; }
  double getDeltaE() const;
  bool getUseMixing() const { return usedmixing_; }

 private:
  Eigen::MatrixXd DensityMatrixGroundState_unres(
      const Eigen::MatrixXd& MOs, Index nocclevels) const;

  void Levelshift(Eigen::MatrixXd& H, const Eigen::MatrixXd& MOs_old,
                  const options& opt, Index nocclevels) const;

  Eigen::MatrixXd BuildErrorMatrix(const Eigen::MatrixXd& dmat,
                                   const Eigen::MatrixXd& H) const;

  double CombinedError(const Eigen::MatrixXd& err_alpha,
                       const Eigen::MatrixXd& err_beta) const;

  options opt_alpha_;
  options opt_beta_;

  AOOverlap* S_ = nullptr;
  Logger* log_ = nullptr;
  Eigen::MatrixXd Sminusahalf;

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
};

}  // namespace xtp
}  // namespace votca

#endif // VOTCA_XTP_UKS_CONVERGENCEACC_H