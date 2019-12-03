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

#pragma once
#ifndef VOTCA_XTP_SYMMETRIC_MATRIX_H
#define VOTCA_XTP_SYMMETRIC_MATRIX_H

#include <iostream>
#include <vector>
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {

/*
 * A symmetric matrix implementation for doubles, acces upper diagonal matrix
 */
class Symmetric_Matrix {
 public:
  Symmetric_Matrix(Index dim) {
    dimension = dim;
    data.resize((dim + 1) * dim / 2);
  }

  Symmetric_Matrix(const Eigen::MatrixXd& full);

  Index size() const { return dimension; }

  double TraceofProd(const Symmetric_Matrix& a) const;

  void AddtoEigenMatrix(Eigen::MatrixXd& full, double factor = 1.0) const;

  void AddtoEigenUpperMatrix(
      Eigen::SelfAdjointView<Eigen::MatrixXd, Eigen::Upper>& upper,
      double factor = 1.0) const;

  Eigen::MatrixXd FullMatrix() const;
  // returns a matrix where only the upper triangle part is filled, rest is set
  // to zero
  Eigen::MatrixXd UpperMatrix() const;

  double& operator()(Index i, Index j) { return data[index(i, j)]; };

  const double& operator()(Index i, Index j) const {
    return data[index(i, j)];
  };

  friend std::ostream& operator<<(std::ostream& out, const Symmetric_Matrix& a);

 private:
  Index index(Index i, Index j) const;

  std::vector<double> data;
  Index dimension;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SYMMETRIC_MATRIX_H
