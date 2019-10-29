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

  _value.clear();

  _grid.clear();

  _coeff.clear();
}

void PadeApprox::addPoint(std::complex<double> frequency, Eigen::MatrixXcd value) {

  this->_grid.push_back(frequency);
  this->_value.push_back(value);
  this->_coeff.push_back(RecursivePolynom(_grid.size() - 1, _grid.size()));
  if (_coeff.at(_coeff.size() - 1).norm() !=
      _coeff.at(_coeff.size() - 1).norm()) {
    std::cout
        << "reject point, unvalid coeficient at w=" << frequency
        << std::endl;  //<<_coeff.at(_coeff.size()-1)<<std::endl<<std::endl;
    _coeff.pop_back();
    _grid.pop_back();
    _value.pop_back();
    _rejected_points++;
  }
}

Eigen::MatrixXcd PadeApprox::RecursivePolynom(int indx, int degree) {

  //Eigen::MatrixXcd temp = Eigen::MatrixXcd(_basis_size, _basis_size);

  if (degree == 1) {
    return _value.at(indx);
  } else {
    Eigen::MatrixXcd temp = RecursivePolynom(indx, degree - 1);

    Eigen::MatrixXcd u = _coeff[degree-2] - temp;
    Eigen::MatrixXcd l = temp * (_grid.at(indx) - _grid.at(degree - 2));
    return u * (l.inverse());
  }
}

Eigen::MatrixXcd PadeApprox::RecursiveA(std::complex<double> frequency, int index) {

  if (index == 0) {
    return Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  } else if (index == 1) {
    return _coeff.at(0);
  } else {
    Eigen::MatrixXcd A = RecursiveA(frequency, index - 1) + (frequency - _grid.at(index - 2)) *
                                                    _coeff.at(index - 1) *
                                                    RecursiveA(frequency, index - 2);
    return A;
  }
}

Eigen::MatrixXcd PadeApprox::RecursiveB(std::complex<double> frequency, int index) {

  if (index == 0) {
    return Eigen::MatrixXcd::Identity(_basis_size, _basis_size);
  } else if (index == 1) {
    return Eigen::MatrixXcd::Identity(_basis_size, _basis_size);
  } else {
    Eigen::MatrixXcd B = RecursiveB(frequency, index - 1) + (frequency - _grid.at(index - 2)) *
                                                    _coeff.at(index - 1) *
                                                    RecursiveB(frequency, index - 2);
    return B;
  }
}

void PadeApprox::initialize(int basis_size) { this->_basis_size = basis_size; }

Eigen::MatrixXcd PadeApprox::evaluatePoint(std::complex<double> frequency) {

  return RecursiveB(frequency, _grid.size()).inverse() *
         RecursiveA(frequency, _grid.size());
}

}  // namespace xtp
}  // namespace votca