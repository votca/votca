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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/cudamatrix.h>

namespace votca {
namespace xtp {
cudaError_t checkCuda(cudaError_t result) {
#if defined(DEBUG)
  if (result != cudaSuccess) {
    std::cerr << "CUDA Runtime Error: " << cudaGetErrorString(result) << "\n";
  }
#endif
  return result;
}

Index count_available_gpus() {
  int count;
  cudaError_t err = cudaGetDeviceCount(&count);
  return 0 ? (err != cudaSuccess) : Index(count);
}

CudaMatrix::CudaMatrix(const Eigen::MatrixXd &matrix,
                       const cudaStream_t &stream)
    : _rows{static_cast<Index>(matrix.rows())},
      _cols{static_cast<Index>(matrix.cols())} {
  _data = alloc_matrix_in_gpu(size_matrix());
  _stream = stream;
  cudaError_t err = cudaMemcpyAsync(_data.get(), matrix.data(), size_matrix(),
                                    cudaMemcpyHostToDevice, stream);
  if (err != 0) {
    throw std::runtime_error("Error copy arrays to device");
  }
}

CudaMatrix::CudaMatrix(Index nrows, Index ncols, const cudaStream_t &stream)
    : _rows{static_cast<Index>(nrows)}, _cols{static_cast<Index>(ncols)} {
  _data = alloc_matrix_in_gpu(size_matrix());
  _stream = stream;
}

CudaMatrix::operator Eigen::MatrixXd() const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(this->rows(), this->cols());
  checkCuda(cudaMemcpyAsync(result.data(), this->data(), this->size_matrix(),
                            cudaMemcpyDeviceToHost, this->_stream));
  return result;
}

void CudaMatrix::copy_to_gpu(const Eigen::MatrixXd &A) {
  size_t size_A = static_cast<Index>(A.size()) * sizeof(double);
  checkCuda(cudaMemcpyAsync(this->data(), A.data(), size_A,
                            cudaMemcpyHostToDevice, _stream));
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