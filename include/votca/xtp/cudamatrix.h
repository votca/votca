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

#include <votca_config.h>
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

inline cudaError_t checkCuda(cudaError_t result) {
#if defined(DEBUG)
  if (result != cudaSuccess) {
    std::cerr << "CUDA Runtime Error: " << cudaGetErrorString(result) << "\n";
  }
#endif
  return result;
}

class CudaMatrix {
 public:
  int size() const { return _rows * _cols; };
  int rows() const { return _rows; };
  int cols() const { return _cols; };
  double *data() const { return _data.get(); };

  CudaMatrix(const Eigen::MatrixXd &matrix, const cudaStream_t &stream)
      : _rows{static_cast<int>(matrix.rows())},
        _cols{static_cast<int>(matrix.cols())} {
    _data = std::move(alloc_matrix_in_gpu(size_matrix()));
    _stream = stream;
    cudaError_t err = cudaMemcpyAsync(_data.get(), matrix.data(), size_matrix(),
                                      cudaMemcpyHostToDevice, stream);
    if (err != 0) {
      throw std::runtime_error("Error copy arrays to device");
    }
  }

  // Allocate memory in the GPU for a matrix
  CudaMatrix(long int nrows, long int ncols)
      : _rows{static_cast<int>(nrows)}, _cols{static_cast<int>(ncols)} {
    _data = std::move(alloc_matrix_in_gpu(size_matrix()));
  }

  // Convert A Cudamatrix to an EigenMatrix
  operator Eigen::MatrixXd() const {
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(this->rows(), this->cols());
    checkCuda(cudaMemcpyAsync(result.data(), this->data(), this->size_matrix(),
                              cudaMemcpyDeviceToHost, this->_stream));
    return result;
  };

  void copy_to_gpu(const Eigen::MatrixXd &A) {
    size_t size_A = static_cast<int>(A.size()) * sizeof(double);
    checkCuda(cudaMemcpyAsync(this->data(), A.data(), size_A,
                              cudaMemcpyHostToDevice, _stream));
  }

  // Unique pointer with custom delete function
  using Unique_ptr_to_GPU_data = std::unique_ptr<double, void (*)(double *)>;

 private:
  Unique_ptr_to_GPU_data alloc_matrix_in_gpu(size_t size_arr) const {
    double *dmatrix;
    throw_if_not_enough_memory_in_gpu(size_arr);
    checkCuda(cudaMalloc(&dmatrix, size_arr));
    Unique_ptr_to_GPU_data dev_ptr(dmatrix,
                                   [](double *x) { checkCuda(cudaFree(x)); });
    return dev_ptr;
  }

  void throw_if_not_enough_memory_in_gpu(size_t requested_memory) const {
    size_t free, total;
    checkCuda(cudaMemGetInfo(&free, &total));

    std::ostringstream oss;
    oss << "There were requested : " << requested_memory
        << "bytes int the device\n";
    oss << "Device Free memory (bytes): " << free
        << "\nDevice total Memory (bytes): " << total << "\n";

    // Raise an error if there is not enough total or free memory in the device
    if (requested_memory > free) {
      oss << "There is not enough memory in the Device!\n";
      throw std::runtime_error(oss.str());
    }
  }

  size_t size_matrix() const { return this->size() * sizeof(double); }

  // Attributes of the matrix in the device
  Unique_ptr_to_GPU_data _data{nullptr,
                               [](double *x) { checkCuda(cudaFree(x)); }};
  cudaStream_t _stream = nullptr;
  int _rows;
  int _cols;
};

}  // namespace xtp
}  // namespace votca

#endif
