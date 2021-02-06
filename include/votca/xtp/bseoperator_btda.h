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

#ifndef VOTCA_XTP_BSEOPERATOR_BTDA_H
#define VOTCA_XTP_BSEOPERATOR_BTDA_H

// Local VOTCA includes
#include "bse_operator.h"
#include "eigen.h"

namespace votca {
namespace xtp {

template <typename MatrixReplacementA, typename MatrixReplacementB>
class HamiltonianOperator;
}
}  // namespace votca
namespace Eigen {
namespace internal {

template <typename MatrixReplacementA, typename MatrixReplacementB>
struct traits<
    votca::xtp::HamiltonianOperator<MatrixReplacementA, MatrixReplacementB>>
    : public Eigen::internal::traits<Eigen::MatrixXd> {};
}  // namespace internal
}  // namespace Eigen

namespace votca {
namespace xtp {
template <typename MatrixReplacementA, typename MatrixReplacementB>
class HamiltonianOperator
    : public Eigen::EigenBase<
          HamiltonianOperator<MatrixReplacementA, MatrixReplacementB>> {
 public:
  // Required typedefs, constants, and method:
  using Scalar = double;
  using RealScalar = double;
  using StorageIndex = votca::Index;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  HamiltonianOperator(const MatrixReplacementA& A, const MatrixReplacementB& B)
      : _A(A), _B(B) {
    _size = 2 * A.cols();
    _diag = get_diagonal();
  };

  Eigen::Index rows() const { return this->_size; }
  Eigen::Index cols() const { return this->_size; }

  template <typename Vtype>
  Eigen::Product<HamiltonianOperator, Vtype, Eigen::AliasFreeProduct> operator*(
      const Eigen::MatrixBase<Vtype>& x) const {
    return Eigen::Product<HamiltonianOperator, Vtype, Eigen::AliasFreeProduct>(
        *this, x.derived());
  }

  Eigen::VectorXd diagonal() const { return _diag; }

  const MatrixReplacementA& _A;
  const MatrixReplacementB& _B;

 private:
  Eigen::VectorXd get_diagonal() const {
    Eigen::VectorXd diag = Eigen::VectorXd::Zero(_size);
    Index half = _size / 2;
    diag.head(half) = _A.diagonal();
    diag.tail(half) = -diag.head(half);
    return diag;
  }

  Index _size;
  Eigen::VectorXd _diag;
};
}  // namespace xtp
}  // namespace votca

namespace Eigen {
namespace internal {

// replacement of the mat*mat operation
template <typename Mtype, typename MatrixReplacementA,
          typename MatrixReplacementB>
struct generic_product_impl<
    votca::xtp::HamiltonianOperator<MatrixReplacementA, MatrixReplacementB>,
    Mtype, DenseShape, DenseShape, GemmProduct>
    : generic_product_impl_base<
          votca::xtp::HamiltonianOperator<MatrixReplacementA,
                                          MatrixReplacementB>,
          Mtype,
          generic_product_impl<votca::xtp::HamiltonianOperator<
                                   MatrixReplacementA, MatrixReplacementB>,
                               Mtype>> {

  typedef typename Product<
      votca::xtp::HamiltonianOperator<MatrixReplacementA, MatrixReplacementB>,
      Mtype>::Scalar Scalar;

  template <typename Dest>
  static void scaleAndAddTo(Dest& dst,
                            const votca::xtp::HamiltonianOperator<
                                MatrixReplacementA, MatrixReplacementB>& op,
                            const Mtype& m, const Scalar& alpha) {
    // returns dst = alpha * op * v
    // alpha must be 1 here
    assert(alpha == Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);

    Index half = op.rows() / 2;
    /**Instead of doing the
    (A   B)*(M1)
    (-B -A) (M2)
     multiplication explicitly for each block
     we reshape M into (M1,M2)
     and multiply A*(M1,M2)
     and then sort the contributions into the resulting vector
     we do the same for B
     * **/
    const Map<const MatrixXd> m_reshaped(m.data(), m.rows() / 2, m.cols() * 2);
    {
      const MatrixXd temp = op._A * m_reshaped;
      const Map<const MatrixXd> temp_unshaped(temp.data(), m.rows(), m.cols());
      dst.topRows(half) = temp_unshaped.topRows(half);
      dst.bottomRows(half) = -temp_unshaped.bottomRows(half);
    }
    {
      const MatrixXd temp = op._B * m_reshaped;
      const Map<const MatrixXd> temp_unshaped2(temp.data(), m.rows(), m.cols());
      dst.topRows(half) += temp_unshaped2.bottomRows(half);
      dst.bottomRows(half) -= temp_unshaped2.topRows(half);
    }
  }
};
}  // namespace internal
}  // namespace Eigen

#endif  // VOTCA_XTP_BSEOPERATOR_BTDA_H
