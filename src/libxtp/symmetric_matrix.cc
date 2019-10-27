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
#include "votca/xtp/symmetric_matrix.h"
#include <iostream>

namespace votca {
namespace xtp {

Symmetric_Matrix::Symmetric_Matrix(const Eigen::MatrixXd& full)
    : dimension(full.rows()) {
  assert(full.rows() == full.cols() && "Input matrix not quadratic");
  data.resize((dimension + 1) * dimension / 2);
  for (long i = 0; i < full.rows(); ++i) {
    for (long j = 0; j <= i; ++j) {
      this->operator()(i, j) = full(i, j);
    }
  }
}

std::ostream& operator<<(std::ostream& out, const Symmetric_Matrix& a) {

  out << "[" << a.dimension << "," << a.dimension << "]\n";
  for (long i = 0; i < a.dimension; ++i) {
    for (long j = 0; j <= i; ++j) {
      out << a(i, j);
      if (i == j) {
        out << "\n";
      } else {
        out << " ";
      }
    }
  }
  return out;
}

double Symmetric_Matrix::TraceofProd(const Symmetric_Matrix& a) const {
  assert(data.size() == a.data.size() && "Matrices do not have same size");
  double result = 0.0;

  for (long i = 0; i < dimension; ++i) {
    const long index = (i * (i + 1)) / 2 + i;
    result += +data[index] * a.data[index];
  }
  for (long i = 0; i < dimension; ++i) {
    const long start = (i * (i + 1)) / 2;
    for (long j = 0; j < i; ++j) {
      result += 2 * data[start + j] * a.data[start + j];
    }
  }
  return result;
}

void Symmetric_Matrix::AddtoEigenMatrix(Eigen::MatrixXd& full,
                                        double factor) const {
  for (long j = 0; j < full.cols(); ++j) {
    const long start = (j * (j + 1)) / 2;
    for (long i = 0; i <= j; ++i) {
      full(i, j) += factor * data[start + i];
    }

    for (long i = j + 1; i < full.rows(); ++i) {
      const long index = (i * (i + 1)) / 2 + j;
      full(i, j) += factor * data[index];
    }
  }
  return;
}

void Symmetric_Matrix::AddtoEigenUpperMatrix(
    Eigen::SelfAdjointView<Eigen::MatrixXd, Eigen::Upper>& upper,
    double factor) const {
  for (long j = 0; j < upper.cols(); ++j) {
    const long start = (j * (j + 1)) / 2;
    for (long i = 0; i <= j; ++i) {
      upper(i, j) += factor * data[start + i];
    }
  }
  return;
}

Eigen::MatrixXd Symmetric_Matrix::FullMatrix() const {
  Eigen::MatrixXd result = Eigen::MatrixXd(dimension, dimension);
  for (long j = 0; j < result.cols(); ++j) {
    const long start = (j * (j + 1)) / 2;
    for (long i = 0; i <= j; ++i) {
      result(i, j) = data[start + i];
    }

    for (long i = j + 1; i < result.rows(); ++i) {
      const long index = (i * (i + 1)) / 2 + j;
      result(i, j) = data[index];
    }
  }
  return result;
}

Eigen::MatrixXd Symmetric_Matrix::UpperMatrix() const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(dimension, dimension);
  for (long j = 0; j < result.cols(); ++j) {
    const long start = (j * (j + 1)) / 2;
    for (long i = 0; i <= j; ++i) {
      result(i, j) = data[start + i];
    }
  }
  return result;
}

long Symmetric_Matrix::Index(long i, long j) const {
  long index;
  if (i >= j) {
    index = (i * (i + 1)) / 2 + j;
  } else {
    index = (j * (j + 1)) / 2 + i;
  }
  return index;
}

}  // namespace xtp
}  // namespace votca
