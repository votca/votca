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
#ifndef VOTCA_XTP_MULTISHIFT_H
#define VOTCA_XTP_MULTISHIFT_H

#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>

namespace votca {
namespace xtp {
class Multishift {
 public:
 
  struct MultiShiftResult{
      std::vector<std::complex<double>> _step_length_a;
      std::vector<std::complex<double>> _step_length_b;
      std::vector<Eigen::VectorXcd> _residue;
      Eigen::VectorXcd _x;
      bool converged=true;    
  };
  
  void setMatrixSize(int size);

  Multishift::MultiShiftResult ComplexBiCG(const Eigen::MatrixXcd& A, const Eigen::VectorXcd& b)const;

  Eigen::VectorXcd DoMultishift(const Eigen::MatrixXcd& A, const Eigen::VectorXcd& b,
                                std::complex<double> w, MultiShiftResult input) const;

 private:
   
  int _matrix_size;

};
}  // namespace xtp
}  // namespace votca
#endif //VOTCA_VOTCA-XTP_MULTISHIFT_H