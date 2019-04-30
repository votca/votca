
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

#include <iostream>
#include <votca/xtp/matrixfreeoperator.h>

namespace votca {
namespace xtp {

Eigen::VectorXd MatrixFreeOperator::diagonal() const {
  Eigen::VectorXd D = Eigen::VectorXd::Zero(_size);
  Eigen::RowVectorXd row_data;
  for (int i = 0; i < _size; i++) {
    row_data = this->row(i);
    D(i) = row_data(i);
  }
  return D;
}

// get the full matrix if we have to
Eigen::MatrixXd MatrixFreeOperator::get_full_matrix() const {
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(_size, _size);
#pragma omp parallel for
  for (int i = 0; i < _size; i++) {
    matrix.row(i) = this->row(i);
  }
  return matrix;
}

// get the size
int MatrixFreeOperator::size() const { return this->_size; }

// set the size
void MatrixFreeOperator::set_size(int size) { this->_size = size; }

}  // namespace xtp
}  // namespace votca