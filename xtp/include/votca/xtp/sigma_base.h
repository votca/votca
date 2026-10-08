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
#ifndef VOTCA_XTP_SIGMA_BASE_H
#define VOTCA_XTP_SIGMA_BASE_H

#include <functional>
#include <string>

// Local VOTCA includes
#include "eigen.h"
#include "screening_kernels.h"

namespace votca {
namespace xtp {

class TCMatrix_gwbse;
class RPA;

class Sigma_base {
 public:
  Sigma_base(TCMatrix_gwbse& Mmn, const RPA& rpa) : Mmn_(Mmn), rpa_(rpa){};

  virtual ~Sigma_base() = default;

  struct options {
    Index homo;
    Index qpmin;
    Index qpmax;
    Index rpamin;
    Index rpamax;
    double eta;
    std::string quadrature_scheme;  // Gaussian-quadrature scheme to use in CDA
    Index order;  // used in numerical integration of CDA Sigma
    double alpha;
  };

  void configure(options opt) {
    opt_ = opt;
    qptotal_ = opt.qpmax - opt.qpmin + 1;
    rpatotal_ = opt.rpamax - opt.rpamin + 1;
  }

  // Calculates full exchange matrix
  Eigen::MatrixXd CalcExchangeMatrix() const;

  /**
   * \brief Static self-energy of an environment reaction field: the
   * COH+SEX pair for v_reac, in the QP window.
   *
   * With the reaction field v_reac = v_12 chi^(2) v_21 of a static
   * environment folded into the screened interaction, W = W_dyn + v_reac,
   * and the frequency-independent part contributes
   *
   *   Sigma^SEX_nn' = - sum_{i occ} (ni|v_reac|n'i)
   *   Sigma^COH_nn' = 1/2 sum_{m}   (nm|v_reac|n'm)
   *
   *   Sigma^reac_nn' = 1/2 sum_m s_m (nm|v_reac|n'm),  s_m = -1 occ, +1 virt,
   *
   * with (nm|v_reac|n'm) = M_nm R M_n'm^T. SEX alone would move the HOMO by
   * twice the polarization energy and leave the LUMO where it is; with COH
   * a level of density psi_n^2 moves by -+ 1/2 (nn|v_reac|nn): occupied
   * levels up, virtual ones down, which is the familiar P_h + P_e gap
   * closing of a polarizable environment.
   *
   * The m-sum runs over the RPA window. COH needs it complete -- it is a
   * closure relation -- so with rpamax below the last MO the result is
   * truncated. v_reac is smooth over the molecule, so that converges fast,
   * but it is not exact.
   *
   * R must be in the FILL-TIME auxiliary frame, R = T^T B T as
   * EnvironmentScreening::SymmetrizedReactionField returns it; it is
   * brought into the frame of the integrals here.
   */
  Eigen::MatrixXd CalcReactionFieldMatrix(const Eigen::MatrixXd& R) const;
  // Calculates correlation diagonal
  Eigen::VectorXd CalcCorrelationDiag(const Eigen::VectorXd& frequencies) const;
  // Calculates correlation off-diagonal (diagonal left zero)
  virtual Eigen::MatrixXd CalcCorrelationOffDiag(
      const Eigen::VectorXd& frequencies) const;

  // Sets up the screening parametrisation
  virtual void PrepareScreening() = 0;

  /// One line about the screening just prepared (empty if nothing to say).
  virtual std::string ScreeningSummary() const { return ""; }
  /// False if CalcCorrelationOffDiag only returns zeros (CDA).
  virtual bool HasOffDiagonal() const { return true; }
  // Calculates Sigma_c diagonal elements
  virtual double CalcCorrelationDiagElementDerivative(
      Index gw_level, double frequency) const = 0;
  virtual double CalcCorrelationDiagElement(Index gw_level,
                                            double frequency) const = 0;
  /// Diagonal element at several frequencies. The default evaluates them one
  /// by one; implementations may do all of them in one pass over the
  /// integrals.
  virtual Eigen::VectorXd CalcCorrelationDiagElements(
      Index gw_level, const Eigen::VectorXd& frequencies) const {
    Eigen::VectorXd result(frequencies.size());
    for (Index i = 0; i < frequencies.size(); ++i) {
      result(i) = CalcCorrelationDiagElement(gw_level, frequencies(i));
    }
    return result;
  }
  // Calculates Sigma_c off-diagonal elements
  virtual double CalcCorrelationOffDiagElement(Index gw_level1, Index gw_level2,
                                               double frequency1,
                                               double frequency2) const = 0;

  /// Memory (bytes) for the level block and the integral chunk of
  /// PairContraction; smaller values only mean more passes.
  void setPairContractionMemory(double block_bytes, double chunk_bytes) {
    block_bytes_ = block_bytes;
    chunk_bytes_ = chunk_bytes;
  }

  /// Times the parts of PrepareScreening in these timings (optional).
  void setTimings(DFTTimings* timings) { timings_ = timings; }

  void ResetDiagEvalCounter() const { diag_eval_counter_.store(0); }
  std::size_t GetDiagEvalCounter() const { return diag_eval_counter_.load(); }

 protected:
  options opt_;
  TCMatrix_gwbse& Mmn_;
  const RPA& rpa_;
  DFTTimings* timings_ = nullptr;

  Index qptotal_ = 0;
  Index rpatotal_ = 0;

  void CountDiagEval() const { diag_eval_counter_.fetch_add(1); }

  /// A(m,n) = sum_{l < nrows, P} X_m(l,P) M_n(l,P) for all pairs of QP
  /// levels m, n, with M_n = Mmn_[n] (QP-window numbering) and X_m from
  /// makeX(m, X) (nrows x auxsize). Evaluated as matrix products over blocks
  /// of m and chunks of the auxiliary index, so the three-centre integrals
  /// are streamed a few times instead of once per pair.
  Eigen::MatrixXd PairContraction(
      Index nrows,
      const std::function<void(Index, Eigen::MatrixXd&)>& makeX) const;

 private:
  double block_bytes_ = 4e9;
  double chunk_bytes_ = 5e8;
  mutable std::atomic<std::size_t> diag_eval_counter_{0};
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SIGMA_BASE_H
