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

#include <complex>
#include <fstream>
#include <math.h>
#include <votca/tools/property.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/multishift.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {

void Multishift::setMatrixSize(int size) {
  this->_matrix_size = size;
}

Multishift::MultiShiftResult Multishift::ComplexBiCG(const Eigen::MatrixXcd& A, const Eigen::VectorXcd& b) const{

  MultiShiftResult result;

  result._residue.push_back(b);
  Eigen::VectorXcd r_t = result._residue[0].conjugate();
  Eigen::VectorXcd p = result._residue[0];
  Eigen::VectorXcd p_t = r_t;

  Eigen::VectorXcd x = Eigen::VectorXcd::Zero(_matrix_size);
 
  double res = 1;
  double tol = 1e-18;
 
  int i = 0;

  int max_iter = 1000000;

  while (res > tol) {
    
    result._step_length_a.push_back(r_t.dot(result._residue[i]) / p_t.dot(A * p));
    x = x + result._step_length_a[i] * p;
    result._residue.push_back(result._residue[i] - result._step_length_a[i] * A * p);
    r_t = r_t - std::conj(result._step_length_a[i]) * A.adjoint() * p_t;
    result._step_length_b.push_back(-1 * (A.adjoint() * p_t).dot(result._residue[i + 1]) / p_t.dot(A * p));
    p = result._residue[i + 1] + result._step_length_b[i] * p;
    p_t = r_t + std::conj(result._step_length_b[i]) * p_t;
    res = result._residue[i].squaredNorm();
    i++;
    if (i == max_iter) {
      
      result.converged=false;
      
      std::cout<<"cBiCG failed"<<std::endl;
      
      x = A.colPivHouseholderQr().solve(b);
    }
  }
  result._x=x;
  return result;
}

Eigen::VectorXcd Multishift::DoMultishift(const Eigen::MatrixXcd& A,
                                          const Eigen::VectorXcd& b,
                                          std::complex<double> w,
                                          MultiShiftResult input) const{

  if (input._residue.empty()) {
    // std::cout<<"Using HouseholderQR"<<std::endl;
    Eigen::MatrixXcd LHS =
        A + w * Eigen::MatrixXcd::Identity(_matrix_size, _matrix_size);
    return LHS.colPivHouseholderQr().solve(b);
  }

  Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(_matrix_size, _matrix_size);

  Eigen::VectorXcd residue = b;
  Eigen::VectorXcd residue_t = residue.conjugate();
  Eigen::VectorXcd direction = residue;
  Eigen::VectorXcd direction_t = residue_t;

  Eigen::VectorXcd x = Eigen::VectorXcd::Zero(_matrix_size);

  std::complex<double> beta;

  std::complex<double> pi = 1;
  std::complex<double> stepsize;
  std::complex<double> pi_p = (1 + w * input._step_length_a[0]) * 1;

  std::complex<double> pi_temp;

  for(int i = 0;i < input._step_length_a.size();i++) {

    stepsize = (pi / pi_p) * input._step_length_a[i];
    x = x + stepsize * direction;
    residue = input._residue[i+ 1] / pi_p;
    residue_t = residue_t - std::conj(input._step_length_a[i]) * (A + w * I).adjoint() * direction_t;

    beta = (pi / pi_p) * (pi / pi_p) * input._step_length_b[i];

    direction = residue + beta * direction;
    direction_t = residue_t + std::conj(beta) * direction_t;

    pi_temp = pi_p;
    pi_p = (1 + w * input._step_length_a[i+1]) * pi_p +
           input._step_length_a[i+1] * input._step_length_b[i+1- 1] / input._step_length_a[i+1- 1] * (pi_p - pi);
    pi = pi_temp;
  }
  
  return x;
}
}  // namespace xtp
}  // namespace votca