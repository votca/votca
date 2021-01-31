/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_MATRIXFREEOPERATOR_H
#define VOTCA_XTP_MATRIXFREEOPERATOR_H

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

class MatrixFreeOperator;
}
}  // namespace votca
namespace Eigen {
namespace internal {
// MatrixReplacement looks-like a Matrix, so let's inherits its traits:
template <>
struct traits<votca::xtp::MatrixFreeOperator>
    : public Eigen::internal::traits<Eigen::MatrixXd> {};
}  // namespace internal
}  // namespace Eigen

namespace votca {
namespace xtp {

class MatrixFreeOperator : public Eigen::EigenBase<MatrixFreeOperator> {
 public:
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  Index rows() const { return this->_size; }
  Index cols() const { return this->_size; }

  template <typename Vtype>
  Eigen::Product<MatrixFreeOperator, Vtype, Eigen::AliasFreeProduct> operator*(
      const Eigen::MatrixBase<Vtype>& x) const {
    return Eigen::Product<MatrixFreeOperator, Vtype, Eigen::AliasFreeProduct>(
        *this, x.derived());
  }

  virtual Eigen::VectorXd diagonal() const = 0;
  virtual Eigen::MatrixXd matmul(const Eigen::MatrixXd& input) const = 0;
  Index size() const;
  void set_size(Index size);

 private:
  Index _size;
};
}  // namespace xtp
}  // namespace votca

namespace Eigen {

namespace internal {

// replacement of the mat*mat operation
template <typename Mtype>
struct generic_product_impl<votca::xtp::MatrixFreeOperator, Mtype, DenseShape,
                            DenseShape, GemmProduct>
    : generic_product_impl_base<
          votca::xtp::MatrixFreeOperator, Mtype,
          generic_product_impl<votca::xtp::MatrixFreeOperator, Mtype>> {

  typedef
      typename Product<votca::xtp::MatrixFreeOperator, Mtype>::Scalar Scalar;

  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const votca::xtp::MatrixFreeOperator& op,
                            const Mtype& m, const Scalar& alpha) {
    // returns dst = alpha * op * v
    // alpha must be 1 here
    assert(alpha == Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);
    dst = op.matmul(m);
  }
};
}  // namespace internal
}  // namespace Eigen

#endif  // VOTCA_XTP_MATRIXFREEOPERATOR_H
