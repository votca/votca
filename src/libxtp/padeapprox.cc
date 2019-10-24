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

#include <fstream>
#include <math.h>
#include <votca/tools/property.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/padeapprox.h>

namespace votca {
namespace xtp {

void PadeApprox::clear() {

  _imval.clear();

  _imgrid.clear();

  _coeff.clear();
}

void PadeApprox::addPoint(std::complex<double> w, Eigen::MatrixXcd val) {

  this->_imgrid.push_back(w);
  this->_imval.push_back(val);
  this->_coeff.push_back(RecursivePolynom(_imgrid.size() - 1, _imgrid.size()));
  if (_coeff.at(_coeff.size() - 1).norm() !=
      _coeff.at(_coeff.size() - 1).norm()) {
    std::cout
        << "reject point, unvalid coeficient at w=" << w
        << std::endl;  //<<_coeff.at(_coeff.size()-1)<<std::endl<<std::endl;
    _coeff.pop_back();
    _imgrid.pop_back();
    _imval.pop_back();
    _rejected_points++;
  }
}

Eigen::MatrixXcd PadeApprox::RecursivePolynom(int indx, int p) {

  //Eigen::MatrixXcd temp = Eigen::MatrixXcd(_basis_size, _basis_size);

  if (p == 1) {
    return _imval.at(indx);
  } else {
    Eigen::MatrixXcd temp = RecursivePolynom(indx, p - 1);

    Eigen::MatrixXcd u = RecursivePolynom(p - 2, p - 1) - temp;
    Eigen::MatrixXcd l = temp * (_imgrid.at(indx) - _imgrid.at(p - 2));
    if (abs(l.determinant()) == 0) {
      std::cout << "det l " << l.determinant() << " p= " << p << std::endl;
      return Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
    }
    return u * (l.inverse());
  }
}

Eigen::MatrixXcd PadeApprox::RecursiveA(std::complex<double> w, int n) {

  if (n == 0) {
    return Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  } else if (n == 1) {
    return _coeff.at(0);
  } else {
    Eigen::MatrixXcd A = RecursiveA(w, n - 1) + (w - _imgrid.at(n - 2)) *
                                                    _coeff.at(n - 1) *
                                                    RecursiveA(w, n - 2);
    return A;
  }
}

Eigen::MatrixXcd PadeApprox::RecursiveB(std::complex<double> w, int n) {

  if (n == 0) {
    return Eigen::MatrixXcd::Identity(_basis_size, _basis_size);
  } else if (n == 1) {
    return Eigen::MatrixXcd::Identity(_basis_size, _basis_size);
  } else {
    Eigen::MatrixXcd B = RecursiveB(w, n - 1) + (w - _imgrid.at(n - 2)) *
                                                    _coeff.at(n - 1) *
                                                    RecursiveB(w, n - 2);
    return B;
  }
}

void PadeApprox::initialize(int basis_size) { this->_basis_size = basis_size; }

Eigen::MatrixXcd PadeApprox::evaluate(std::complex<double> w) {

  if (abs(RecursiveB(w, _imgrid.size()).determinant()) == 0) {
    std::cout << std::endl
              << "Warning, Matrix close to singular, det="
              << abs(RecursiveB(w, _imgrid.size()).determinant()) << std::endl;
    return Eigen::MatrixXcd(_basis_size, _basis_size);
  }
  return RecursiveB(w, _imgrid.size()).inverse() *
         RecursiveA(w, _imgrid.size());
}

void PadeApprox::test() {

  std::cout << "Started Pade Test" << std::endl;

  initialize(20);

  // define Grid
  std::vector<std::complex<double>> grid;

  std::complex<double> i(0, 1);

  grid.push_back(1);
  grid.push_back(2);
  grid.push_back(3);
  grid.push_back(4);
  grid.push_back(5);
  grid.push_back(6);
  grid.push_back(7);
  grid.push_back(8);

  // std::cout<<"Gridpoints are = "<<grid<<std::endl;
  for (int i = 0; i < grid.size(); i++) {
    std::cout << "Gridpoint: " << grid.at(i) << std::endl;
  }
  std::vector<Eigen::MatrixXcd> val;

  // val.push_back(Eigen::MatrixXcd::Ones(2,2)+Eigen::MatrixXcd::Identity(2,2));
  // val.push_back(Eigen::MatrixXcd::Ones(2,2));
  //        val.push_back(Eigen::MatrixXcd::Ones(2,2)+3*Eigen::MatrixXcd::Identity(2,2));
  //        val.push_back(Eigen::MatrixXcd::Ones(2,2)+4*Eigen::MatrixXcd::Identity(2,2));
  //        val.push_back(Eigen::MatrixXcd::Ones(2,2)+5*Eigen::MatrixXcd::Identity(2,2));

  val.push_back(0.000001 * Eigen::MatrixXcd::Random(20, 20));
  val.push_back(0.000001 * Eigen::MatrixXcd::Random(20, 20));
  val.push_back(0.000001 * Eigen::MatrixXcd::Random(20, 20));
  val.push_back(0.000001 * Eigen::MatrixXcd::Random(20, 20));
  val.push_back(0.000001 * Eigen::MatrixXcd::Random(20, 20));
  val.push_back(0.000001 * Eigen::MatrixXcd::Random(20, 20));
  val.push_back(0.000001 * Eigen::MatrixXcd::Random(20, 20));
  val.push_back(0.000001 * Eigen::MatrixXcd::Random(20, 20));
  val.push_back(0.000001 * Eigen::MatrixXcd::Random(20, 20));
  // val.push_back(0.0001*Eigen::MatrixXcd::Random(20,20));
  for (int i = 0; i < grid.size(); i++) {
    std::cout << "Value: " << val.at(i)(1, 1) << std::endl;
  }
  for (int i = 0; i < grid.size(); i++) {

    addPoint(grid.at(i), val.at(i));
  }
  //        for(int i=0; i<_coeff.size(); i++){
  //
  //            std::cout<<"i= "<<i<<"coeff size "<<_coeff.size()<<std::endl;
  //
  //            std::cout<<"Printing Coefficients "<<_coeff.at(i)<<std::endl;
  //
  //        }

  std::cout << "Done calculating coefficients" << std::endl;
  //        for(int i=0; i<grid.size()-_rejected_points; i++){
  //            std::cout<<"Coeff: "<<i<<std::endl<<_coeff.at(i)<<std::endl;
  //        }

  // std::cout<<"f(0)="<<evaluate(0)(1,1)<<std::endl;
  //        std::cout<<"f(0.25)="<<evaluate(0.25)(1,1)<<std::endl;
  //        std::cout<<"f(0.5)="<<evaluate(0.5)(1,1)<<std::endl;
  //        std::cout<<"f(0.75)="<<evaluate(0.75)(1,1)<<std::endl;
  std::cout << "f(1)=" << evaluate(1)(1, 1) << std::endl;
  //        std::cout<<"f(1.25)="<<evaluate(1.25)(1,1)<<std::endl;
  //        std::cout<<"f(1.5)="<<evaluate(1.5)(1,1)<<std::endl;
  //        std::cout<<"f(1.75)="<<evaluate(1.75)(1,1)<<std::endl;
  std::cout << "f(2)=" << evaluate(1.5)(1, 1) << std::endl;
  //        std::cout<<"f(2.25)="<<evaluate(2.25)(1,1)<<std::endl;
  //        std::cout<<"f(2.5)="<<evaluate(2.5)(1,1)<<std::endl;
  //        std::cout<<"f(2.75)="<<evaluate(2.75)(1,1)<<std::endl;
  std::cout << "f(3)=" << evaluate(3)(1, 1) << std::endl;
  //        std::cout<<"f(3.25)="<<evaluate(3.25)(1,1)<<std::endl;
  //        std::cout<<"f(3.5)="<<evaluate(3.5)(1,1)<<std::endl;
  //        std::cout<<"f(3.75)="<<evaluate(3.75)(1,1)<<std::endl;
  std::cout << "f(4)=" << evaluate(4)(1, 1) << std::endl;
  //        std::cout<<"f(4.25)="<<evaluate(4.25)(1,1)<<std::endl;
  //        std::cout<<"f(4.5)="<<evaluate(4.5)(1,1)<<std::endl;
  //        std::cout<<"f(4.75)="<<evaluate(4.75)(1,1)<<std::endl;
  std::cout << "f(5)=" << evaluate(5)(1, 1) << std::endl;
}
}  // namespace xtp
}  // namespace votca