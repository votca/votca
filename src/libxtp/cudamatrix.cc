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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/cudamatrix.h"

namespace votca {
namespace xtp {
void checkCuda(cudaError_t result) {
  if (result != cudaSuccess) {
    throw std::runtime_error(std::string("CUDA Runtime Error: ") +
                             cudaGetErrorString(result));
  }
}

void checkCublas(cublasStatus_t result) {
  if (result != CUBLAS_STATUS_SUCCESS) {
    throw std::runtime_error(std::string("CUBLAS Runtime Error: ") +
                             cudaGetErrorEnum(result));
  }
}

std::string cudaGetErrorEnum(cublasStatus_t error) {
  switch (error) {
    case CUBLAS_STATUS_SUCCESS:
      return "CUBLAS_STATUS_SUCCESS";
    case CUBLAS_STATUS_NOT_INITIALIZED:
      return "CUBLAS_STATUS_NOT_INITIALIZED";
    case CUBLAS_STATUS_ALLOC_FAILED:
      return "CUBLAS_STATUS_ALLOC_FAILED";
    case CUBLAS_STATUS_INVALID_VALUE:
      return "CUBLAS_STATUS_INVALID_VALUE";
    case CUBLAS_STATUS_ARCH_MISMATCH:
      return "CUBLAS_STATUS_ARCH_MISMATCH";
    case CUBLAS_STATUS_MAPPING_ERROR:
      return "CUBLAS_STATUS_MAPPING_ERROR";
    case CUBLAS_STATUS_EXECUTION_FAILED:
      return "CUBLAS_STATUS_EXECUTION_FAILED";
    case CUBLAS_STATUS_INTERNAL_ERROR:
      return "CUBLAS_STATUS_INTERNAL_ERROR";
    case CUBLAS_STATUS_NOT_SUPPORTED:
      return "CUBLAS_STATUS_NOT_SUPPORTED";
    case CUBLAS_STATUS_LICENSE_ERROR:
      return "CUBLAS_STATUS_LICENSE_ERROR";
  }
  return "<unknown>";
}

Index count_available_gpus() {
  int count;
  cudaError_t err = cudaGetDeviceCount(&count);
  return (err != cudaSuccess) ? 0 : Index(count);
}

CudaMatrix::CudaMatrix(Index nrows, Index ncols, const cudaStream_t &stream)
    : _ld(nrows), _cols(ncols) {
  _data = alloc_matrix_in_gpu(size_matrix());
  _stream = stream;
}

CudaMatrix::operator Eigen::MatrixXd() const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(this->rows(), this->cols());
  checkCuda(cudaMemcpyAsync(result.data(), this->data(), this->size_matrix(),
                            cudaMemcpyDeviceToHost, this->_stream));
  checkCuda(cudaStreamSynchronize(this->_stream));
  return result;
}

CudaMatrix::Unique_ptr_to_GPU_data CudaMatrix::alloc_matrix_in_gpu(
    size_t size_arr) const {
  double *dmatrix;
  throw_if_not_enough_memory_in_gpu(size_arr);
  checkCuda(cudaMalloc(&dmatrix, size_arr));
  Unique_ptr_to_GPU_data dev_ptr(dmatrix,
                                 [](double *x) { checkCuda(cudaFree(x)); });
  return dev_ptr;
}

void CudaMatrix::throw_if_not_enough_memory_in_gpu(
    size_t requested_memory) const {
  size_t free, total;
  checkCuda(cudaMemGetInfo(&free, &total));

  std::ostringstream oss;
  oss << "There were requested : " << requested_memory
      << "bytes Index the device\n";
  oss << "Device Free memory (bytes): " << free
      << "\nDevice total Memory (bytes): " << total << "\n";

  // Raise an error if there is not enough total or free memory in the device
  if (requested_memory > free) {
    oss << "There is not enough memory in the Device!\n";
    throw std::runtime_error(oss.str());
  }
}

}  // namespace xtp
}  // namespace votca
