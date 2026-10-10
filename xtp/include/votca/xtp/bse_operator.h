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

// Standard includes
#include <memory>

// Local VOTCA includes
#include "eigen.h"
#include "matrixfreeoperator.h"
#include "symmetric_tiles.h"
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

  /// Keep the screened direct term (Hd for cd, Hd2 for cd2; unit
  /// prefactor) in this matrix, which is built on first use unless it
  /// already holds it. Operators of the same BSE can share it (e.g. the
  /// singlet and triplet A, or Hd for the analysis); null disables it.
  void set_direct_cache(std::shared_ptr<SymmetricTiles> cache) {
    direct_cache_ = std::move(cache);
  }
  bool direct_term_cached() const {
    return direct_cache_ != nullptr && direct_cache_->built();
  }
  /// Tag of the direct term an operator caches: 1 for Hd, 2 for Hd2
  static constexpr Index DirectKind() {
    return (cd != 0) ? 1 : ((cd2 != 0) ? 2 : 0);
  }
  /// Without the dense cache the screened direct term is applied in large
  /// products, one block of Hamiltonian rows (all v1 for one c1) at a time;
  /// true selects the row-by-row kernel (one Hamiltonian row at a time).
  void use_row_kernel(bool rows) { row_kernel_ = rows; }

 private:
  // the screened direct term (cd or cd2 part, unit prefactor) as the lower
  // triangle of a symmetric matrix
  void BuildDirectMatrix() const;
  // the screened direct term applied block by block (see use_row_kernel)
  Eigen::MatrixXd ApplyDirectBlocks(const Eigen::MatrixXd& input) const;
  Eigen::MatrixXd ApplyDirectRows(const Eigen::MatrixXd& input) const;

  std::shared_ptr<SymmetricTiles> direct_cache_;
  bool row_kernel_ = false;
  // stacked slices for the block kernel: for cd the vv blocks of the
  // occupied slices (row v1*n_v + v2), for cd2 the cv blocks of the virtual
  // slices (row c1*n_v + v2)
  mutable Eigen::MatrixXd stacked_;
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
