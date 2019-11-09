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
#ifndef VOTCA_XTP_EIGEN_H
#define VOTCA_XTP_EIGEN_H

#include <votca/tools/eigen.h>
#include <votca/tools/types.h>
#include <votca/xtp/votca_config.h>
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Array<votca::Index, Eigen::Dynamic, 1> ArrayXl;

namespace votca {
namespace xtp {

// Stores matrix and energy together
class Mat_p_Energy {
 public:
  Mat_p_Energy(Index rows, Index cols)
      : _energy(0.0), _matrix(Eigen::MatrixXd::Zero(rows, cols)){};
  Mat_p_Energy(double e, const Eigen::MatrixXd& mat)
      : _energy(e), _matrix(mat){};
  Mat_p_Energy(double e, Eigen::MatrixXd&& mat)
      : _energy(e), _matrix(std::move(mat)){};

  Index rows() const { return _matrix.rows(); }
  Index cols() const { return _matrix.cols(); }
  Eigen::MatrixXd& matrix() { return _matrix; }
  double& energy() { return _energy; }
  const Eigen::MatrixXd& matrix() const { return _matrix; }
  double energy() const { return _energy; }

 private:
  double _energy;
  Eigen::MatrixXd _matrix;
};

namespace OPENMP {
inline Index getMaxThreads() {
  Index nthreads = 1;
#ifdef _OPENMP
  nthreads = Index(omp_get_max_threads());
#endif
  return nthreads;
}

inline Index getThreadId() {
  Index thread_id = 0;
#ifdef _OPENMP
  thread_id = Index(omp_get_thread_num());
#endif
  return thread_id;
}

#ifdef _OPENMP
inline void setMaxThreads(Index threads) {
  if (threads > 0) {
    omp_set_num_threads(int(threads));
  }
}
#else
inline void setMaxThreads(Index) {}
#endif
}  // namespace OPENMP
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EIGEN_H
