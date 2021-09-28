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
#ifndef VOTCA_XTP_INCREMENTALFOCKBUILDER_H
#define VOTCA_XTP_INCREMENTALFOCKBUILDER_H

#include "votca/xtp/logger.h"
#include <votca/tools/types.h>
namespace votca {
namespace xtp {

// Small Wrapper class to build incremental fock matrix
class IncrementalFockBuilder {
 public:
  IncrementalFockBuilder(Logger& log, double start_threshold,
                         Index fock_matrix_reset)
      : log_(log),
        start_incremental_F_threshold_(start_threshold),
        fock_matrix_reset_(fock_matrix_reset) {}

  void Configure(const Eigen::MatrixXd& dmat) {
    Ddiff_ = dmat;
    Dlast_ = dmat;
  }

  void Start(Index iteration, double DiisError) {
    if (!incremental_Fbuild_started_ &&
        DiisError < start_incremental_F_threshold_) {
      incremental_Fbuild_started_ = true;
      reset_incremental_fock_formation_ = false;
      last_reset_iteration_ = iteration - 1;
      next_reset_threshold_ = DiisError / 10.0;
      XTP_LOG(Log::error, log_)
          << TimeStamp() << " Using incremental 4c build from here"
          << std::flush;
    }
  }

  void resetMatrices(Eigen::MatrixXd& J, Eigen::MatrixXd& K,
                     const Eigen::MatrixXd& dmat) {
    if (reset_incremental_fock_formation_ || !incremental_Fbuild_started_) {
      J.setZero();
      K.setZero();
      Ddiff_ = dmat;
    }
  }

  const Eigen::MatrixXd& getDmat_diff() const { return Ddiff_; }

  void UpdateCriteria(double DiisError, Index Iteration) {
    if (reset_incremental_fock_formation_ && incremental_Fbuild_started_) {
      reset_incremental_fock_formation_ = false;
      last_reset_iteration_ = Iteration;
      next_reset_threshold_ = DiisError / 10.0;
      XTP_LOG(Log::error, log_)
          << TimeStamp() << " Reset incremental 4c build" << std::flush;
    }
  }

  void UpdateDmats(const Eigen::MatrixXd& dmat, double DiisError,
                   Index Iteration) {
    if (DiisError < next_reset_threshold_ ||
        Iteration - last_reset_iteration_ > fock_matrix_reset_) {
      reset_incremental_fock_formation_ = true;
    }
    Ddiff_ = dmat - Dlast_;
    Dlast_ = dmat;
  }

 private:
  Logger& log_;
  double start_incremental_F_threshold_;  // Diis error from which to start
                                          // using incremental builds
  Index fock_matrix_reset_;  // After how many iterations the fock matrix should
                             // be reset regardless

  Eigen::MatrixXd Ddiff_;
  Eigen::MatrixXd Dlast_;

  bool reset_incremental_fock_formation_ = false;
  bool incremental_Fbuild_started_ = false;
  double next_reset_threshold_ = 0.0;
  Index last_reset_iteration_ = 0;
};  // namespace xtp

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_INCREMENTALFOCKBUILDER_H
