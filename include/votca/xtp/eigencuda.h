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

#include <cublas_v2.h>
#include <curand.h>
#include <sstream>
#include <vector>
#include <votca/xtp/eigen.h>

/*
 * \brief Perform Tensor-matrix multiplications in a GPU
 *
 * The `EigenCuda` class handles the allocation and deallocation of arrays on
 * the GPU.
 */

namespace votca {
namespace xtp {

inline cudaError_t checkCuda(cudaError_t result) {
// Check Cuda error
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    std::cerr << "CUDA Runtime Error: " << cudaGetErrorString(result) << "\n";
  }
#endif
  return result;
}

// Structure with the sizes to call ?GEMM
struct ShapesOfMatrices {
  int A_rows;
  int A_cols;
  int B_rows;
  int B_cols;
  int C_rows;

  ShapesOfMatrices(long int a_rows, long int a_cols, long int b_rows,
                   long int b_cols, long int crows)
      : A_rows{static_cast<int>(a_rows)},
        A_cols{static_cast<int>(a_cols)},
        B_rows{static_cast<int>(b_rows)},
        B_cols{static_cast<int>(b_cols)},
        C_rows{static_cast<int>(crows)} {}
};

// col Major for CUDA
using Mat =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

class EigenCuda {
 public:
  EigenCuda() {
    cublasCreate(&_handle);
    cudaStreamCreate(&_stream);
  }
  ~EigenCuda();

  EigenCuda(const EigenCuda &) = delete;
  EigenCuda &operator=(const EigenCuda &) = delete;

  // Perform the triple matrix multiplication A * matrix * C, for the vector
  // of matrices given by tensor
  std::vector<Mat> matrix_tensor_matrix_mult(const Mat &A,
                                             const std::vector<Mat> &tensor,
                                             const Mat &C);

  // Perform a multiplication between a matrix and a tensor
  std::vector<Mat> right_matrix_tensor_mult(const std::vector<Mat> &tensor,
                                            const Mat &A) const;

 private:
  void check_available_memory_in_gpu(size_t required) const;
  void alloc_mem_in_gpu(double **x, std::size_t n) const;
  void free_mem_in_gpu(double *x) const;
  double *alloc_matrix_in_gpu(size_t size_matrix) const;

  // Allocate memory for a matrix and copy it to the device
  double *copy_matrix_to_gpu(const Mat &matrix) const;

  // Invoke the ?gemm function of cublas
  void gemm(ShapesOfMatrices shapes, const double *dA, const double *dB,
            double *dC) const;

  // The cublas handles allocates hardware resources on the host and device.
  cublasHandle_t _handle;

  // Asynchronous stream
  cudaStream_t _stream;
};

}  // namespace xtp
}  // namespace votca

#endif /*XTP_EIGEN_CUDA__H */
