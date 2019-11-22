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

#ifndef __XTP_CUDA_MATRIX__H
#define __XTP_CUDA_MATRIX__H

#include <votca/xtp/votca_config.h>
#ifndef USE_CUDA
#error Cuda not enabled
#endif

#include <cublas_v2.h>
#include <curand.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <votca/xtp/eigen.h>

/*
 * \brief Matrix Representation inside an Nvidia GPU
 */
namespace votca {
namespace xtp {

cudaError_t checkCuda(cudaError_t result);

Index count_available_gpus();

class CudaMatrix {
 public:
  Index size() const { return _rows * _cols; };
  Index rows() const { return _rows; };
  Index cols() const { return _cols; };
  double *data() const { return _data.get(); };

  CudaMatrix(const Eigen::MatrixXd &matrix, const cudaStream_t &stream);

  // Allocate memory in the GPU for a matrix
  CudaMatrix(Index nrows, Index ncols, const cudaStream_t &stream);

  // Convert A Cudamatrix to an EigenMatrix
  operator Eigen::MatrixXd() const;

  void copy_to_gpu(const Eigen::MatrixXd &A);

 private:
  // Unique pointer with custom delete function
  using Unique_ptr_to_GPU_data = std::unique_ptr<double, void (*)(double *)>;

  Unique_ptr_to_GPU_data alloc_matrix_in_gpu(size_t size_arr) const;

  void throw_if_not_enough_memory_in_gpu(size_t requested_memory) const;

  size_t size_matrix() const { return this->size() * sizeof(double); }

  // Attributes of the matrix in the device
  Unique_ptr_to_GPU_data _data{nullptr,
                               [](double *x) { checkCuda(cudaFree(x)); }};
  cudaStream_t _stream = nullptr;
  Index _rows;
  Index _cols;
};

}  // namespace xtp
}  // namespace votca

#endif
