/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef __XTP_EIGEN_CUDA__H
#define __XTP_EIGEN_CUDA__H
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cublas_v2.h>
#include <curand.h>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace votca {
namespace xtp {

// Structure with the sizes to call ?GEMM
struct Shapes {
  int A_rows;
  int A_cols;
  int B_rows;
  int B_cols;
  int C_rows;

  Shapes(long int _a_rows, long int _a_cols, long int _b_rows, long int _b_cols,
         long int _c_rows)
      : A_rows{static_cast<int>(_a_rows)},
        A_cols{static_cast<int>(_a_cols)},
        B_rows{static_cast<int>(_b_rows)},
        B_cols{static_cast<int>(_b_cols)},
        C_rows{static_cast<int>(_c_rows)} {}
};
// col Major for CUDA
template <typename T>
using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

template <typename T>
class EigenCuda {

 public:
  EigenCuda() { cublasCreate(&_handle); }
  EigenCuda(bool pinned) : _pinned{pinned} { cublasCreate(&_handle); }

  // Deallocate both the handler and allocated arrays
  ~EigenCuda();

  // Remove the copy operations
  EigenCuda(const EigenCuda &) = delete;
  EigenCuda &operator=(const EigenCuda &) = delete;

  // Matrix matrix multiplication
  Mat<T> dot(Mat<T> &A, Mat<T> &B);

  // Perform the triple matrix multiplication A * matrix * C, for the vector
  // of matrices given by tensor
  std::vector<Mat<T>> triple_tensor_product(Mat<T> &A, Mat<T> &C,
                                            std::vector<Mat<T>> &tensor);

 private:
  // Allocate memory in the device
  void fun_alloc(T **x, std::size_t n) const;

  // Deallocate memory from the device
  void fun_free(T *x) const;

  // Copy matricex to the device
  unsigned initialize_Matrix(Mat<T> &A, bool copy_to_device = true);

  // Invoke the ?gemm function of cublas
  void gemm(Shapes shapes, std::tuple<unsigned, unsigned, unsigned> ids);

  // Deallocate certain matrix from the device
  void free_matrix(unsigned id);

  // Cuda variables
  cublasHandle_t _handle;
  bool _pinned = false;

  // Allocation booking
  unsigned _counter = 0;
  std::unordered_map<unsigned, T *> _allocated;

  // Scalar constanst for calling blas
  T _alpha = 1.;
  T _beta = 0.;
  const T *_pa = &_alpha;
  const T *_pb = &_beta;
};

}  // namespace xtp
}  // namespace votca

#endif /*XTP_EIGEN_CUDA__H */
