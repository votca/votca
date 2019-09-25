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
// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
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
      
    
  class InnerIterator {
   public:
    InnerIterator(const HamiltonianOperator& xpr, const Eigen::Index& id)
        : _xpr(xpr), _id(id){};

    InnerIterator& operator++() {
      _row++;
      return *this;
    }
    operator bool() const {
      return _row < _xpr._size;
    }  // DO not use the size method, it returns linear dimension*linear
       // dimension i.e _size^2
    double value() const { return _xpr(_row, _id); }
    Eigen::Index row() const { return _row; }
    Eigen::Index col() const { return _id; }
    Eigen::Index index() const { return row(); }

   private:
    const HamiltonianOperator& _xpr;
    const Eigen::Index _id;
    Eigen::Index _row = 0;
  };

  Eigen::Index rows() const { return this->_size; }
  Eigen::Index cols() const { return this->_size; }
  Eigen::Index outerSize() const { return this->_size; }

  template <typename Vtype>
  Eigen::Product<HamiltonianOperator, Vtype, Eigen::AliasFreeProduct>
      operator*(const Eigen::MatrixBase<Vtype>& x) const {
    return Eigen::Product<HamiltonianOperator, Vtype,
                          Eigen::AliasFreeProduct>(*this, x.derived());
  }

  // this is not a fast method
  const double& operator()(const size_t i, const size_t j) const {
    if (i==j) {
      return _diag(i);
    } else {
      Eigen::RowVectorXd row_out = row(i);
      return  row_out(j);
    }
  };

  //  get a row of the operator
  Eigen::RowVectorXd row(int index) const {
    /* Returna row of the operator 
    H = [ A   B
         -B* -A* ]
    */
    int lsize = this->_size;
    int halfsize = lsize/2;
   
    Eigen::RowVectorXd row_out = Eigen::RowVectorXd::Zero(lsize);

    if (index < halfsize)
    {

      Eigen::RowVectorXd a = this->_A.row(index);
      Eigen::RowVectorXd b = this->_B.row(index);    
      row_out.head(halfsize) = a;
      row_out.tail(halfsize) = b;

    } else {

      Eigen::RowVectorXd a = this->_A.row(index-halfsize);
      Eigen::RowVectorXd b = this->_B.row(index-halfsize);    
      row_out.head(halfsize) = -b;
      row_out.tail(halfsize) = -a;
    }
    return row_out;
  }

  Eigen::VectorXd get_diagonal() const {
    Eigen::VectorXd diag = Eigen::VectorXd::Zero(_size);
#pragma omp parallel for
    for (int i=0; i < _size; i++) {
      Eigen::RowVectorXd r = row(i);
      diag(i) = r(i);
    }
    return diag;
  }

  // get the full matrix if we have to
  Eigen::MatrixXd get_full_matrix() const {
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(_size, _size);
  
    for (int i = 0; i < _size; i++) {
      matrix.row(i) = row(i);
    }
    return matrix;
  }

 private:

  const MatrixReplacementA & _A;
  const MatrixReplacementB & _B;
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

// make the mat vect product
#pragma omp parallel for
    for (int i = 0; i < op.rows(); i++) {
      dst(i) = op.row(i) * v;
    }
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

// make the mat mat product
#pragma omp parallel for
    for (int i = 0; i < op.rows(); i++) {
      const Eigen::RowVectorXd row = op.row(i) * m;
      dst.row(i) = row;
    }
  }
};
}
}

#endif  //__VOTCA_BSEOP_BTDA_H