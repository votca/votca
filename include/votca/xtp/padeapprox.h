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

#pragma once
#ifndef PADEAPPROX_H
#define PADEAPPROX_H

#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>

namespace votca {
namespace xtp {
class PadeApprox {
 public:
  PadeApprox(){};

  Eigen::MatrixXcd evaluate(std::complex<double> w);

  void addPoint(std::complex<double> w, Eigen::MatrixXcd val);

  void initialize(int basis_size);

  void clear();

  void test();

 private:
  std::vector<Eigen::MatrixXcd> _imval;

  std::vector<std::complex<double>> _imgrid;

  std::vector<Eigen::MatrixXcd> _coeff;

  int _basis_size;

  int _rejected_points = 0;

  Eigen::MatrixXcd RecursivePolynom(int indx, int p);

  Eigen::MatrixXcd RecursiveA(std::complex<double> w, int n);

  Eigen::MatrixXcd RecursiveB(std::complex<double> w, int n);
};
}  // namespace xtp
}  // namespace votca
#endif