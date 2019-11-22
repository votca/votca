
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

Eigen::RowVectorXd MatrixFreeOperator::OperatorRow(Index) const {
  return Eigen::RowVectorXd::Zero(0);
}

Eigen::MatrixXd MatrixFreeOperator::OperatorBlock(long, long) const {
  return Eigen::MatrixXd::Zero(0, 0);
}

Eigen::VectorXd MatrixFreeOperator::diagonal() const {
  Eigen::VectorXd D = Eigen::VectorXd::Zero(_size);
  if (useRow()) {
#pragma omp parallel for schedule(guided)
    for (Index i = 0; i < _size; i++) {
      Eigen::RowVectorXd row_data = this->OperatorRow(i);
      D(i) = row_data(i);
    }
  }

  if (useBlock()) {
    Index blocksize = getBlocksize();
    if (size() % blocksize != 0) {
      throw std::runtime_error("blocksize is not a multiple of matrix size");
    }
    Index blocks = size() / blocksize;

// this is inefficient if blocks<num_ofthreads
#pragma omp parallel for schedule(guided)
    for (Index i = 0; i < blocks; i++) {
      Eigen::MatrixXd block = OperatorBlock(i, i);
      D.segment(i * blocksize, blocksize) += block.diagonal();
    }
  }

  return D;
}

// get the full matrix if we have to
Eigen::MatrixXd MatrixFreeOperator::get_full_matrix() const {
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(_size, _size);
  if (useRow()) {
#pragma omp parallel for schedule(guided)
    for (Index i = 0; i < _size; i++) {
      matrix.row(i) = this->OperatorRow(i);
    }
  }

  if (useBlock()) {
    Index blocksize = getBlocksize();
    if (size() % blocksize != 0) {
      throw std::runtime_error("blocksize is not a multiple of matrix size");
    }
    Index blocks = size() / blocksize;

// this is inefficient if blocks<num_ofthreads
#pragma omp parallel for schedule(guided)
    for (Index i_row = 0; i_row < blocks; i_row++) {
      for (Index i_col = 0; i_col < blocks; i_col++) {
        matrix.block(i_row * blocksize, i_col * blocksize, blocksize,
                     blocksize) += OperatorBlock(i_row, i_col);
      }
    }
  }

  return matrix;
}

// get the size
long MatrixFreeOperator::size() const { return this->_size; }

// set the size
void MatrixFreeOperator::set_size(Index size) { this->_size = size; }

}  // namespace xtp
}  // namespace votca
