/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#ifndef _VOTCA_XTP_BSE_OPERATOR_H
#define _VOTCA_XTP_BSE_OPERATOR_H

#include <votca/xtp/eigen.h>
#include <votca/xtp/matrixfreeoperator.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

struct BSEOperator_Options {
  int homo;
  int rpamin;
  int qpmin;
  int vmin;
  int cmax;
};

template <int cqp, int cx, int cd, int cd2>
class BSE_OPERATOR : public MatrixFreeOperator {

 public:
  BSE_OPERATOR(const Eigen::VectorXd& Hd_operator, const TCMatrix_gwbse& Mmn,
               const Eigen::MatrixXd& Hqp)
      : _epsilon_0_inv(Hd_operator), _Mmn(Mmn), _Hqp(Hqp){};

  void configure(BSEOperator_Options opt);

  Eigen::RowVectorXd OperatorRow(int index) const override;

  bool useBlock() const override { return cx != 0; }

  int getBlocksize() const override { return _bse_ctotal; }

  Eigen::MatrixXd OperatorBlock(int row, int col) const override;

 private:
  Eigen::RowVectorXd Hqp_row(int index) const;
  Eigen::RowVectorXd Hd_row(int index) const;
  Eigen::RowVectorXd Hd2_row(int index) const;
  Eigen::MatrixXd HxBlock(int row, int col) const;

  BSEOperator_Options _opt;
  int _bse_size;
  int _bse_vtotal;
  int _bse_ctotal;
  int _bse_cmin;

  const Eigen::VectorXd& _epsilon_0_inv;
  const TCMatrix_gwbse& _Mmn;
  const Eigen::MatrixXd& _Hqp;
};

// type defs for the different operators
typedef BSE_OPERATOR<1, 2, 1, 0> SingletOperator_TDA;
typedef BSE_OPERATOR<1, 0, 1, 0> TripletOperator_TDA;

typedef BSE_OPERATOR<1, 4, 1, 1> SingletOperator_BTDA_ApB;
typedef BSE_OPERATOR<1, 0, 1, 1> TripletOperator_BTDA_ApB;
typedef BSE_OPERATOR<1, 0, 1, -1> Operator_BTDA_AmB;

typedef BSE_OPERATOR<0, 2, 0, 1> SingletOperator_BTDA_B;

typedef BSE_OPERATOR<1, 0, 0, 0> HqpOperator;
typedef BSE_OPERATOR<0, 1, 0, 0> HxOperator;
typedef BSE_OPERATOR<0, 0, 1, 0> HdOperator;
typedef BSE_OPERATOR<0, 0, 0, 1> Hd2Operator;

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_BSE_OP_H */
