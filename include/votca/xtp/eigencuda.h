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

#ifndef __XTP_EIGEN_CUDA__H
#define __XTP_EIGEN_CUDA__H

#include <cublas_v2.h>
#include <curand.h>
#include <iostream>
#include <memory>
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
#if defined(DEBUG)
  if (result != cudaSuccess) {
    std::cerr << "CUDA Runtime Error: " << cudaGetErrorString(result) << "\n";
  }
#endif
  return result;
}

// Data of the matrix stored in the GPU
class CudaMatrix {
 public:
  int size() const { return _rows * _cols; };
  int rows() const { return _rows; };
  int cols() const { return _cols; };
  double *pointer() const { return _pointer.get(); };

  CudaMatrix(long int nrows, long int ncols)
      : _rows{static_cast<int>(nrows)}, _cols{static_cast<int>(ncols)} {
    size_t size_matrix = _rows * _cols * sizeof(double);
    _pointer = std::move(alloc_matrix_in_gpu(size_matrix));
  }

  // Allocate memory for a matrix and copy it to the device
  void copy_matrix_to_gpu(const Eigen::MatrixXd &matrix,
                          const cudaStream_t &_stream) const;

  // Unique pointer with custom delete function
  using double_unique_ptr = std::unique_ptr<double, void (*)(double *)>;

 private:
  double_unique_ptr alloc_matrix_in_gpu(size_t size_matrix) const;

  double_unique_ptr _pointer{nullptr,
                             [](double *x) { checkCuda(cudaFree(x)); }};
  int _rows;
  int _cols;
};

void free_mem_in_gpu(double *x);

/* \brief The EigenCuda class offload Eigen operations to an *Nvidia* GPU using
 * the CUDA language.
 * The Cublas handle is the context manager for all the resources needed by
 * Cublas. While a stream is a queue of sequential operations executed in the
 * Nvidia device.
 */
class EigenCuda {
 public:
  EigenCuda() {
    cublasCreate(&_handle);
    cudaStreamCreate(&_stream);
  }
  ~EigenCuda();

  EigenCuda(const EigenCuda &) = delete;
  EigenCuda &operator=(const EigenCuda &) = delete;

  // Perform a multiplication between a matrix and a tensor
  void right_matrix_tensor_mult(std::vector<Eigen::MatrixXd> &tensor,
                                const Eigen::MatrixXd &A) const;

  // Perform matrix1 * matrix2 * matrix3 multiplication
  Eigen::MatrixXd triple_matrix_mult(const CudaMatrix &A,
                                     const Eigen::MatrixXd &matrix,
                                     const CudaMatrix &C) const;

  const cudaStream_t &get_stream() const { return _stream; };

 private:
  void throw_if_not_enough_memory_in_gpu(size_t required) const;

  // Invoke the ?gemm function of cublas
  void gemm(const CudaMatrix &A, const CudaMatrix &B, CudaMatrix &C) const;

  // The cublas handles allocates hardware resources on the host and device.
  cublasHandle_t _handle;

  // Asynchronous stream
  cudaStream_t _stream;
};

}  // namespace xtp
}  // namespace votca

#endif /*XTP_EIGEN_CUDA__H */
