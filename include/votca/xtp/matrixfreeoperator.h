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

#pragma once
#ifndef __VOTCA_TOOLS_MATRIX_FREE_OPERATOR_H
#define __VOTCA_TOOLS_MATRIX_FREE_OPERATOR_H
#include <votca/xtp/eigen.h>

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

  // convenience function
  Eigen::MatrixXd get_full_matrix() const;
  Eigen::VectorXd diagonal() const;
  long size() const;
  void set_size(long size);

  virtual bool useRow() const { return true; }
  virtual bool useBlock() const { return false; }

  virtual long getBlocksize() const { return 0; }

  // extract row/col of the operator
  virtual Eigen::RowVectorXd OperatorRow(long index) const;

  virtual Eigen::MatrixXd OperatorBlock(long row, long col) const;

 private:
  long _size;
};
}  // namespace xtp
}  // namespace votca

namespace Eigen {

namespace internal {

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
    if (op.useRow()) {
// make the mat vect product
#pragma omp parallel for schedule(guided)
      for (long i = 0; i < op.rows(); i++) {
        dst(i) = op.OperatorRow(i) * v;
      }
    }

    if (op.useBlock()) {
      long blocksize = op.getBlocksize();
      if (op.size() % blocksize != 0) {
        throw std::runtime_error("blocksize is not a multiple of matrix size");
      }
      long blocks = op.size() / blocksize;

// this is inefficient if blocks<num_ofthreads
#pragma omp parallel for schedule(guided)
      for (long i_row = 0; i_row < blocks; i_row++) {
        for (long i_col = 0; i_col < blocks; i_col++) {
          dst.segment(i_row * blocksize, blocksize) +=
              op.OperatorBlock(i_row, i_col) *
              v.segment(i_col * blocksize, blocksize);
        }
      }
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
    if (op.useRow()) {
#pragma omp parallel for
      for (long i = 0; i < op.rows(); i++) {
        const Eigen::RowVectorXd row = op.OperatorRow(i) * m;
        dst.row(i) = row;
      }
    }

    if (op.useBlock()) {
      long blocksize = op.getBlocksize();
      if (op.size() % blocksize != 0) {
        throw std::runtime_error("blocksize is not a multiple of matrix size");
      }
      long blocks = op.size() / blocksize;
      // this uses the fact that all our matrices are symmetric, i.e we can
      // reuse half the blocks
      std::vector<Eigen::MatrixXd> thread_wiseresult(
          votca::xtp::OPENMP::getMaxThreads(),
          Eigen::MatrixXd::Zero(dst.rows(), dst.cols()));
#pragma omp parallel for schedule(guided)
      for (long i_row = 0; i_row < blocks; i_row++) {
        for (long i_col = i_row; i_col < blocks; i_col++) {
          Eigen::MatrixXd blockmat = op.OperatorBlock(i_row, i_col);
          thread_wiseresult[votca::xtp::OPENMP::getThreadId()].block(
              i_row * blocksize, 0, blocksize, dst.cols()) +=
              blockmat * m.block(i_col * blocksize, 0, blocksize, m.cols());
          if (i_row != i_col) {
            thread_wiseresult[votca::xtp::OPENMP::getThreadId()].block(
                i_col * blocksize, 0, blocksize, dst.cols()) +=
                blockmat.transpose() *
                m.block(i_row * blocksize, 0, blocksize, m.cols());
          }
        }
      }
      for (const Eigen::MatrixXd mat : thread_wiseresult) {
        dst += mat;
      }
    }
  }
};
}  // namespace internal
}  // namespace Eigen

#endif  //__VOTCA_TOOLS_MATRIX_FREE_OPERATOR_H
