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
      : epsilon_0_inv_(Hd_operator), Mmn_(Mmn), Hqp_(Hqp){};

  void configure(BSEOperator_Options opt);

  // This method sets up the diagonal of the hermitian BSE hamiltonian.
  // Otherwise see the matmul function
  Eigen::VectorXd diagonal() const;
  /*
   * This is the main routine for setting up the hermitian parts of the BSE
   * For the non-hermitian case look at bseoperator_btda.h which combines
   * the bse_operators. The operator is never set up explicitly, instead only
   * the product of it with an input matrix is computed. using the template
   * arguements Different parts of the hamiltonian can be constructed. In
   * general it is inefficient to set them up independently (unless you need it
   * for analysis) thus the function combines all parts.
   */
  Eigen::MatrixXd matmul(const Eigen::MatrixXd& input) const;

  /// Keep the screened direct term as a dense matrix, built on first use, if
  /// it needs at most this many bytes (0, the default, disables it).
  void set_direct_cache_limit(double bytes) { direct_cache_limit_ = bytes; }
  bool direct_term_cached() const { return direct_built_; }
  /// Without the dense cache the screened direct term is applied in large
  /// products, one block of Hamiltonian rows (all v1 for one c1) at a time;
  /// true selects the older row-by-row kernel (always used with GPUs).
  void use_row_kernel(bool rows) { row_kernel_ = rows; }

 private:
  // the screened direct term (cd or cd2 part, with prefactor) as a dense
  // matrix, row-major so each Hamiltonian row is contiguous
  void BuildDirectMatrix() const;
  // the screened direct term applied block by block (see use_row_kernel)
  Eigen::MatrixXd ApplyDirectBlocks(const Eigen::MatrixXd& input) const;
  Eigen::MatrixXd ApplyDirectRows(const Eigen::MatrixXd& input) const;

  double direct_cache_limit_ = 0.0;
  bool row_kernel_ = false;
  mutable bool direct_built_ = false;
  // stacked slices for the block kernel: for cd the vv blocks of the
  // occupied slices (row v1*n_v + v2), for cd2 the cv blocks of the virtual
  // slices (row c1*n_v + v2)
  mutable Eigen::MatrixXd stacked_;
  mutable Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      direct_;
  BSEOperator_Options opt_;
  Index bse_size_;
  Index bse_vtotal_;
  Index bse_ctotal_;
  Index bse_cmin_;

  const Eigen::VectorXd& epsilon_0_inv_;
  const TCMatrix_gwbse& Mmn_;
  const Eigen::MatrixXd& Hqp_;
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
