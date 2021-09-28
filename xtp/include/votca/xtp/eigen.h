/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// CMake Generated file
// clang-format off
// order seems to matter here
#include "votca_xtp_config.h"
#include <votca/tools/votca_tools_config.h>
//clang-format on

// VOTCA includes
#include <votca/tools/eigen.h>
#include <votca/tools/types.h>
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Array<votca::Index, Eigen::Dynamic, 1> ArrayXl;

namespace votca {
namespace xtp {

inline bool XTP_HAS_MKL_OVERLOAD() {

  bool mkl_overload = false;
#ifdef EIGEN_USE_MKL_ALL
  mkl_overload = true;
#endif
  bool mkl_found = false;
#ifdef MKL_FOUND
  mkl_found = true;
#endif
  if (mkl_overload && mkl_found) {
    return true;
  } else {
    return false;
  }
}

// Stores matrix and energy together
class Mat_p_Energy {
 public:
    Mat_p_Energy()
      :  energy_(0.0),  matrix_(Eigen::MatrixXd::Zero(0, 0)){};

  Mat_p_Energy(Index rows, Index cols)
      :  energy_(0.0),  matrix_(Eigen::MatrixXd::Zero(rows, cols)){};
  Mat_p_Energy(double e, const Eigen::MatrixXd& mat)
      :  energy_(e),  matrix_(mat){};
  Mat_p_Energy(double e, Eigen::MatrixXd&& mat)
      :  energy_(e),  matrix_(std::move(mat)){};

  Mat_p_Energy operator+(const Mat_p_Energy& other) const {
    Mat_p_Energy result = *this;
    result. energy_ += other. energy_;
    result. matrix_ += other. matrix_;
    return result;
  }

  Index rows() const { return  matrix_.rows(); }
  Index cols() const { return  matrix_.cols(); }
  Eigen::MatrixXd& matrix() { return  matrix_; }
  double& energy() { return  energy_; }
  const Eigen::MatrixXd& matrix() const { return  matrix_; }
  double energy() const { return  energy_; }

 private:
  double  energy_;
  Eigen::MatrixXd  matrix_;
};

//Stores the diadicProduct of a vector with itself
class AxA {
   public:
    AxA(const Eigen::Vector3d& a) {
       data_.segment<3>(0) = a.x() * a;
       data_.segment<2>(3) = a.y() * a.segment<2>(1);
       data_[5] = a.z() * a.z();
    }
    inline const double& xx() const { return  data_[0]; }
    inline const double& xy() const { return  data_[1]; }
    inline const double& xz() const { return  data_[2]; }
    inline const double& yy() const { return  data_[3]; }
    inline const double& yz() const { return  data_[4]; }
    inline const double& zz() const { return  data_[5]; }

   private:
    Eigen::Matrix<double, 6, 1>  data_;
  };




#pragma omp declare reduction (+:Mat_p_Energy: omp_out=omp_out+omp_in)\
     initializer(omp_priv=Mat_p_Energy(omp_orig.rows(),omp_orig.cols()))

#pragma omp declare reduction (+: Eigen::VectorXd: omp_out=omp_out+omp_in)\
     initializer(omp_priv=Eigen::VectorXd::Zero(omp_orig.size()))

#pragma omp declare reduction (+: Eigen::MatrixXd: omp_out=omp_out+omp_in)\
     initializer(omp_priv=Eigen::MatrixXd::Zero(omp_orig.rows(),omp_orig.cols()))

#pragma omp declare reduction (+: Eigen::Matrix3d: omp_out=omp_out+omp_in)\
     initializer(omp_priv=Eigen::Matrix3d::Zero())

#pragma omp declare reduction (+: Eigen::Vector3d: omp_out=omp_out+omp_in)\
     initializer(omp_priv=Eigen::Vector3d::Zero())

namespace OPENMP {
inline Index getMaxThreads() {
  Index nthreads = 1;
#ifdef  _OPENMP
  nthreads = Index(omp_get_max_threads());
#endif
  return nthreads;
}

inline bool InsideActiveParallelRegion(){
#ifdef  _OPENMP
  return omp_in_parallel();
#endif
return false;
}

inline Index getThreadId() {
  Index thread_id = 0;
#ifdef  _OPENMP
  thread_id = Index(omp_get_thread_num());
#endif
  return thread_id;
}

#ifdef  _OPENMP
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

#endif // VOTCA_XTP_EIGEN_H
