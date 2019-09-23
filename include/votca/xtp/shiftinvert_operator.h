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
#ifndef _VOTCA_XTP_SHIFTINVERT_OPERATOR_H
#define _VOTCA_XTP_SHIFTINVERT_OPERATOR_H

#include <votca/xtp/eigen.h>
#include <votca/xtp/matrixfreeoperator.h>
#include <votca/xtp/threecenter.h>
#include <Eigen/IterativeLinearSolvers>


namespace votca {
  namespace xtp {

template <typename MatrixReplacement>
class ShiftInvertOperator : public MatrixFreeOperator {

public:

  // constructor
  ShiftInvertOperator(const MatrixReplacement &A) :
    _Aop(A), _size(A.cols()) {};

  // get size 
  Eigen::Index rows() const {return this->_size;}
  Eigen::Index cols() const {return this->_size;}

  //  get a row of the operator
  Eigen::RowVectorXd row(int index) const {
    return this->_Aop.row(index);
  }

  // set the shift that is alwyas null in our case ..
  void set_shift(double sigma)
  {
    linear_solver.compute(_Aop);
  }

  // shift invert operation
  void perform_op(const double *x_in, double *y_out)
  {
    Eigen::Map<const Eigen::VectorXd> x(x_in,_size);
    Eigen::Map<Eigen::VectorXd> y(y_out,_size);
    y.noalias() = linear_solver.solve(x);
  }

private:

  MatrixReplacement _Aop;
  double _lambda;
  int _size;
  Eigen::BiCGSTAB<MatrixReplacement> linear_solver;

};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_SSHIFTINVERT_OPERATOR_H */
