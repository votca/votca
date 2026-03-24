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

#pragma once
#ifndef VOTCA_XTP_BSE_OPERATOR_UKS_H
#define VOTCA_XTP_BSE_OPERATOR_UKS_H

#include "eigen.h"
#include "matrixfreeoperator.h"
#include "threecenter.h"

namespace votca {
namespace xtp {

struct BSEOperatorUKS_Options {
  Index homo_alpha;
  Index homo_beta;
  Index rpamin;
  Index qpmin;
  Index vmin;
  Index cmax;
};

template <Index cqp, Index cx, Index cd, Index cd2>
class BSE_OPERATOR_UKS final : public MatrixFreeOperator {
 public:
  BSE_OPERATOR_UKS(const Eigen::VectorXd& epsilon_0_inv,
                   const TCMatrix_gwbse_spin& Mmn,
                   const Eigen::MatrixXd& Hqp_alpha,
                   const Eigen::MatrixXd& Hqp_beta)
      : epsilon_0_inv_(epsilon_0_inv),
        Mmn_(Mmn),
        Hqp_alpha_(Hqp_alpha),
        Hqp_beta_(Hqp_beta) {}

  void configure(BSEOperatorUKS_Options opt);

  Eigen::VectorXd diagonal() const override;
  Eigen::MatrixXd matmul(const Eigen::MatrixXd& input) const override;
  Eigen::MatrixXd dense_matrix() const;

 private:
  struct SpinBlockInfo {
    Index homo = 0;
    Index vmin_rpa = 0;
    Index cmin_rpa = 0;
    Index vtotal = 0;
    Index ctotal = 0;
    Index size = 0;
    Index offset = 0;
  };

  BSEOperatorUKS_Options opt_;
  SpinBlockInfo alpha_;
  SpinBlockInfo beta_;

  Index size_total_ = 0;

  const Eigen::VectorXd& epsilon_0_inv_;
  const TCMatrix_gwbse_spin& Mmn_;
  const Eigen::MatrixXd& Hqp_alpha_;
  const Eigen::MatrixXd& Hqp_beta_;

  void setup_block(SpinBlockInfo& blk, Index homo, Index offset);

  Eigen::VectorXd Hqp_row(const Eigen::MatrixXd& Hqp, const SpinBlockInfo& blk,
                          Index v1, Index c1) const;

  void add_qp_block(Eigen::MatrixXd& y, const Eigen::MatrixXd& x,
                    const SpinBlockInfo& blk, const Eigen::MatrixXd& Hqp) const;

  void add_exchange_block(Eigen::MatrixXd& y, const Eigen::MatrixXd& x,
                          const SpinBlockInfo& out_blk,
                          const SpinBlockInfo& in_blk,
                          const TCMatrix_gwbse& Mout, const TCMatrix_gwbse& Min,
                          double prefactor) const;

  void add_direct_block(Eigen::MatrixXd& y, const Eigen::MatrixXd& x,
                        const SpinBlockInfo& out_blk,
                        const SpinBlockInfo& in_blk, const TCMatrix_gwbse& Mout,
                        const TCMatrix_gwbse& Min, double prefactor) const;

  void add_direct2_block(Eigen::MatrixXd& y, const Eigen::MatrixXd& x,
                         const SpinBlockInfo& out_blk,
                         const SpinBlockInfo& in_blk,
                         const TCMatrix_gwbse& Mout, const TCMatrix_gwbse& Min,
                         double prefactor) const;

  void add_direct_cross_tda_block(Eigen::MatrixXd& y, const Eigen::MatrixXd& x,
                                  const SpinBlockInfo& out_blk,
                                  const SpinBlockInfo& in_blk,
                                  const TCMatrix_gwbse& Mout,
                                  const TCMatrix_gwbse& Min,
                                  double prefactor) const;
};

// TDA A block: Hqp + Hx - Hd
typedef BSE_OPERATOR_UKS<1, 1, 1, 0> ExcitonUKSOperator_TDA;

// full BSE B block: Hx - Hd2
typedef BSE_OPERATOR_UKS<0, 1, 0, 1> ExcitonUKSOperator_BTDA_B;

// direct-only operators used for perturbative dynamical screening
typedef BSE_OPERATOR_UKS<0, 0, 1, 0> HdUKSOperator;
typedef BSE_OPERATOR_UKS<0, 0, 0, 1> Hd2UKSOperator;

typedef BSE_OPERATOR_UKS<1, 0, 0, 0> QpUKSOperator;
typedef BSE_OPERATOR_UKS<0, 1, 0, 0> HxUKSOperator;

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BSE_OPERATOR_UKS_H