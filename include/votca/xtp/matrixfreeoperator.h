/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_TOOLS_MATRIX_FREE_OPERATOR_H
#define __VOTCA_TOOLS_MATRIX_FREE_OPERATOR_H
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {

class MatrixFreeOperator : public Eigen::EigenBase<Eigen::MatrixXd> {
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

  // convenience function
  Eigen::MatrixXd get_full_matrix() const;
  Eigen::VectorXd diagonal() const;
  int size() const;
  void set_size(int size);

  // extract row/col of the operator
  virtual Eigen::VectorXd col(int index) const = 0;

 private:
  int _size;
  Eigen::VectorXd diag_el;

 private:
};
}  // namespace xtp
}  // namespace votca

// user defined reduction
#pragma omp declare reduction(+: Eigen::VectorXd: omp_out=omp_out+omp_in )\
  initializer(omp_priv=Eigen::VectorXd::Zero(omp_orig.size()))

namespace Eigen {

namespace internal {

template <>
struct traits<votca::xtp::MatrixFreeOperator>
    : public Eigen::internal::traits<Eigen::MatrixXd> {};

// replacement of the mat*vect operation
template <typename Vtype>
struct generic_product_impl<votca::xtp::MatrixFreeOperator, Vtype, DenseShape,
                            DenseShape, GemvProduct>
    : generic_product_impl_base<
          votca::xtp::MatrixFreeOperator, Vtype,
          generic_product_impl<votca::xtp::MatrixFreeOperator, Vtype>> {

  typedef
      typename Product<votca::xtp::MatrixFreeOperator, Vtype>::Scalar Scalar;

  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const votca::xtp::MatrixFreeOperator& op,
                            const Vtype& v, const Scalar& alpha) {
    // returns dst = alpha * op * v
    // alpha must be 1 here
    assert(alpha == Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);

// make the mat vect product
#pragma omp parallel for reduction(+ : dst)
    for (int i = 0; i < op.cols(); i++) {
      dst += v(i) * op.col(i);
    }
  }
};

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

// make the mat mat product
#pragma omp parallel for
    for (int j = 0; j < m.cols(); j++) {
      for (int i = 0; i < op.cols(); i++) dst.col(j) += m(i, j) * op.col(i);
    }
  }
};
}  // namespace internal
}  // namespace Eigen

#endif  //__VOTCA_TOOLS_MATRIX_FREE_OPERATOR_H