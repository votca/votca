///*
// *            Copyright 2009-2019 The VOTCA Development Team
// *                       (http://www.votca.org)
// *
// *      Licensed under the Apache License, Version 2.0 (the "License")
// *
// * You may not use this file except in compliance with the License.
// * You may obtain a copy of the License at
// *
// *              http://www.apache.org/licenses/LICENSE-2.0
// *
// * Unless required by applicable law or agreed to in writing, software
// * distributed under the License is distributed on an "AS IS" BASIS,
// * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// * See the License for the specific language governing permissions and
// * limitations under the License.
// */

///*
//* Implementation of the Pade approximation for matrix valued functions,
//* according to the algorithm proposed in the paper "Solving the Eliashberg 
//* equations by means of n-point Pade Approximation" by H.J. Vidberg and J.W. 
//* Serene (1997).
//*/


#pragma once
#ifndef VOTCA_XTP_PADEAPPROX_H
#define VOTCA_XTP_PADEAPPROX_H

#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>

namespace votca {
namespace xtp {
class PadeApprox {
 public:
  PadeApprox(){};

  Eigen::MatrixXcd evaluatePoint(std::complex<double> frequency);

  void addPoint(std::complex<double> frequency, Eigen::MatrixXcd value);

  void initialize(int basis_size);

  void clear();

 private:
  std::vector<Eigen::MatrixXcd> _value;

  std::vector<std::complex<double>> _grid;

  std::vector<Eigen::MatrixXcd> _coeff;

  int _basis_size;

  int _rejected_points = 0;

  Eigen::MatrixXcd RecursivePolynom(int indx, int degree);

  Eigen::MatrixXcd RecursiveA(std::complex<double> frequency, int index);

  Eigen::MatrixXcd RecursiveB(std::complex<double> frequency, int index);
};
}  // namespace xtp
}  // namespace votca
#endif