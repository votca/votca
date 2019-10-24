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

void Multishift::setBasisSize(double basis_size) {
  this->_basis_size = basis_size;
}

Eigen::VectorXcd Multishift::CBiCG(Eigen::MatrixXcd A, Eigen::VectorXcd b) {

  _r.clear();
  _a.clear();
  _b.clear();

  _r.push_back(b);
  Eigen::VectorXcd r_t = _r.at(0).conjugate();
  Eigen::VectorXcd p = _r.at(0);
  Eigen::VectorXcd p_t = r_t;

  Eigen::VectorXcd x = Eigen::VectorXcd::Zero(_basis_size);

  double res = 1;
  double tol = 0.0000001;

  int i = 0;

  int max_iter = 10000;

  while (res > tol) {

    _a.push_back(r_t.dot(_r.at(i)) / p_t.dot(A * p));
    x = x + _a.at(i) * p;
    _r.push_back(_r.at(i) - _a.at(i) * A * p);
    r_t = r_t - std::conj(_a.at(i)) * A.adjoint() * p_t;
    _b.push_back(-1 * (A.adjoint() * p_t).dot(_r.at(i + 1)) / p_t.dot(A * p));
    p = _r.at(i + 1) + _b.at(i) * p;
    p_t = r_t + std::conj(_b.at(i)) * p_t;
    // res=(A*x-b).norm();
    res = _r.at(i).squaredNorm();
    i++;
    if (i == max_iter) {
      // std::cout<<"Max iter reached, res="<<res<<std::endl<<"Using
      // HousholderQR instead."<<std::endl;

      _r.clear();
      _a.clear();
      _b.clear();

      return A.colPivHouseholderQr().solve(b);
    }
  }
  return x;
}

Eigen::VectorXcd Multishift::DoMultishift(Eigen::MatrixXcd A,
                                          Eigen::VectorXcd b,
                                          std::complex<double> w) {

  if (_r.empty()) {
    // std::cout<<"Using HouseholderQR"<<std::endl;
    Eigen::MatrixXcd LHS =
        A + w * Eigen::MatrixXcd::Identity(_basis_size, _basis_size);
    return LHS.colPivHouseholderQr().solve(b);
  }

  Eigen::MatrixXcd I;
  I.setIdentity(_basis_size, _basis_size);

  Eigen::VectorXcd r = b;
  Eigen::VectorXcd r_t = r.conjugate();
  Eigen::VectorXcd p = r;
  Eigen::VectorXcd p_t = r_t;

  Eigen::VectorXcd x = Eigen::VectorXcd::Zero(_basis_size);

  std::complex<double> beta;

  std::complex<double> pi = 1;
  std::complex<double> a;
  std::complex<double> pi_p = (1 + w * _a.at(0)) * 1;

  std::complex<double> pi_temp;

  int i = 0;

  while (i < _a.size()) {

    a = (pi / pi_p) * _a.at(i);
    x = x + a * p;
    r = _r.at(i + 1) / pi_p;
    r_t = r_t - std::conj(_a.at(i)) * (A + w * I).adjoint() * p_t;

    beta = (pi / pi_p) * (pi / pi_p) * _b.at(i);

    p = r + beta * p;
    p_t = r_t + std::conj(beta) * p_t;

    i++;

    if (i >= _a.size()) {
      break;
    }
    pi_temp = pi_p;
    pi_p = (1 + w * _a.at(i)) * pi_p +
           _a.at(i) * _b.at(i - 1) / _a.at(i - 1) * (pi_p - pi);
    pi = pi_temp;
  }
  return x;
}

void Multishift::testMultishift() {

  Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(_basis_size, _basis_size);

  Eigen::MatrixXcd I;

  I.setIdentity(_basis_size, _basis_size);

  std::cout << "A=" << std::endl << A << std::endl;

  Eigen::VectorXcd b = Eigen::VectorXcd::Random(_basis_size);

  std::cout << "b=" << std::endl << b << std::endl;

  Eigen::VectorXcd x = CBiCG(A, b);

  std::cout << "x=" << std::endl << x << std::endl;

  Eigen::VectorXcd res = (A)*x - b;

  std::cout << "res=" << std::endl << res << std::endl;

  if (res.norm() < 0.001) {
    std::cout << "cBiCG Test successful" << std::endl;
  }

  std::complex<double> w(1, 0);

  Eigen::VectorXcd x_w = DoMultishift(A, b, w);

  std::cout << "x_w=" << std::endl << x_w << std::endl;

  Eigen::VectorXcd res_w = (A + w * I) * x_w - b;

  std::cout << "res_w=" << std::endl << res_w << std::endl;

  if (res_w.norm() < 0.001) {
    std::cout << "Multishift Test successful" << std::endl;
  }
}

}  // namespace xtp
}  // namespace votca