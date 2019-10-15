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


#ifndef __VOTCA_BSEOP_BTDA_H
#define __VOTCA_BSEOP_BTDA_H
#include <votca/xtp/bse_operator.h>
#include <votca/xtp/eigen.h>
namespace votca {
namespace xtp {

template<typename MatrixReplacementA, typename MatrixReplacementB>
class HamiltonianOperator;

}
}  // namespace votca
namespace Eigen {
namespace internal {

template <typename MatrixReplacementA, typename MatrixReplacementB>
struct traits<votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>>
    : public Eigen::internal::traits<Eigen::MatrixXd> {};
}  // namespace internal
}  // namespace Eigen

namespace votca {
namespace xtp {
template <typename MatrixReplacementA, typename MatrixReplacementB>
class HamiltonianOperator
    : public Eigen::EigenBase<HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>> {
 public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  HamiltonianOperator(const MatrixReplacementA &A,
                      const MatrixReplacementB &B) 
  : _A(A), _B(B) 
  { 
    _size = 2*A.cols() ;
    _diag = get_diagonal();
   };

  Eigen::Index rows() const { return this->_size; }
  Eigen::Index cols() const { return this->_size; }
  
  template <typename Vtype>
  Eigen::Product<HamiltonianOperator, Vtype, Eigen::AliasFreeProduct>
      operator*(const Eigen::MatrixBase<Vtype>& x) const {
    return Eigen::Product<HamiltonianOperator, Vtype,
                          Eigen::AliasFreeProduct>(*this, x.derived());
  }

  Eigen::VectorXd get_diagonal() const {
    Eigen::VectorXd diag = Eigen::VectorXd::Zero(_size);
    assert(_size%2==0);
    int half = _size/2;
    diag.head(half) = _A.diagonal();
    diag.tail(half) = -diag.head(half);
    return diag;
  }

  Eigen::VectorXd diagonal() const {return _diag;}
  
  // get the full matrix if we have to
  Eigen::MatrixXd get_full_matrix() const {
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(_size, _size);
    int half = _size/2;
    matrix.topLeftCorner(half,half) = _A.get_full_matrix();
    matrix.topRightCorner(half,half) = _B.get_full_matrix();
    matrix.bottomLeftCorner(half,half) = -matrix.topRightCorner(half,half);
    matrix.bottomRightCorner(half,half) = -matrix.topLeftCorner(half,half);
    return matrix;
  }


  const MatrixReplacementA & _A;
  const MatrixReplacementB & _B;

 private:
  int _size;
  Eigen::VectorXd _diag; 

};
}  // namespace xtp
}  // namespace votca

namespace Eigen {
namespace internal {

// replacement of the mat*vect operation
template <typename Vtype, typename MatrixReplacementA, typename MatrixReplacementB>
struct generic_product_impl<votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>, Vtype, DenseShape,
                            DenseShape, GemvProduct>
    : generic_product_impl_base<
          votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>, Vtype,
          generic_product_impl<votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>, Vtype> >{

  typedef
      typename Product<votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>, Vtype>::Scalar Scalar;

  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>& op,
                            const Vtype& v, const Scalar& alpha) {
    // returns dst = alpha * op * v
    // alpha must be 1 here
    assert(alpha == Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);

    int half = op.rows()/2;
    dst.head(half) = op._A*v.head(half) + op._B*v.tail(half);
    dst.tail(half) = - (op._B*v.head(half) + op._A*v.tail(half));

  }
};

// replacement of the mat*mat operation
template <typename Mtype, typename MatrixReplacementA, typename MatrixReplacementB>
struct generic_product_impl<votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>, 
                            Mtype, DenseShape, DenseShape, GemmProduct>
    : generic_product_impl_base<
          votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>, Mtype,
          generic_product_impl<votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>, Mtype>> {

  typedef
      typename Product<votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>, Mtype>::Scalar Scalar;

  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const votca::xtp::HamiltonianOperator<MatrixReplacementA,MatrixReplacementB>& op,
                            const Mtype& m, const Scalar& alpha) {
    // returns dst = alpha * op * v
    // alpha must be 1 here
    assert(alpha == Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);

    int half = op.rows()/2;
    dst.topRows(half) = op._A * m.topRows(half) + op._B * m.bottomRows(half);
    dst.bottomRows(half) = - (op._B * m.topRows(half) + (op._A * m.bottomRows(half)));

  }
};
}
}

#endif  //__VOTCA_BSEOP_BTDA_H