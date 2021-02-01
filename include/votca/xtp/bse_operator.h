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
#ifndef VOTCA_XTP_BSE_OPERATOR_H
#define VOTCA_XTP_BSE_OPERATOR_H

// Local VOTCA includes
#include "eigen.h"
#include "matrixfreeoperator.h"
#include "threecenter.h"

namespace votca {
namespace xtp {

struct BSEOperator_Options {
  Index homo;
  Index rpamin;
  Index qpmin;
  Index vmin;
  Index cmax;
};

template <Index cqp, Index cx, Index cd, Index cd2>
class BSE_OPERATOR final : public MatrixFreeOperator {

 public:
  BSE_OPERATOR(const Eigen::VectorXd& Hd_operator, const TCMatrix_gwbse& Mmn,
               const Eigen::MatrixXd& Hqp)
      : _epsilon_0_inv(Hd_operator), _Mmn(Mmn), _Hqp(Hqp){};

  void configure(BSEOperator_Options opt);

  Eigen::VectorXd diagonal() const;
  Eigen::MatrixXd matmul(const Eigen::MatrixXd& input) const;

 private:
  Eigen::RowVectorXd Hqp_row(Index v1, Index c1) const;
  Eigen::MatrixXd HxBlock(Index v1, Index v2) const;

  BSEOperator_Options _opt;
  Index _bse_size;
  Index _bse_vtotal;
  Index _bse_ctotal;
  Index _bse_cmin;

  const Eigen::VectorXd& _epsilon_0_inv;
  const TCMatrix_gwbse& _Mmn;
  const Eigen::MatrixXd& _Hqp;
};

// type defs for the different operators
typedef BSE_OPERATOR<1, 2, 1, 0> SingletOperator_TDA;
typedef BSE_OPERATOR<1, 0, 1, 0> TripletOperator_TDA;

typedef BSE_OPERATOR<0, 2, 0, 1> SingletOperator_BTDA_B;

typedef BSE_OPERATOR<1, 0, 0, 0> HqpOperator;
typedef BSE_OPERATOR<0, 1, 0, 0> HxOperator;
typedef BSE_OPERATOR<0, 0, 1, 0> HdOperator;
typedef BSE_OPERATOR<0, 0, 0, 1> Hd2Operator;

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BSE_OPERATOR_H
