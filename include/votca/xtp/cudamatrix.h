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

#include <cublas_v2.h>
#include <curand.h>
#include <memory>
#include <votca/xtp/eigen.h>

/*
 * \brief Matrix Representation inside an Nvidia GPU
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

class CudaMatrix {
 public:
  int size() const { return _rows * _cols; };
  int rows() const { return _rows; };
  int cols() const { return _cols; };
  double *pointer() const { return _pointer.get(); };

  CudaMatrix(const Eigen::MatrixXd &matrix, const cudaStream_t &stream)
      : _rows{static_cast<int>(matrix.rows())},
        _cols{static_cast<int>(matrix.cols())} {
    size_t size_matrix = this->size() * sizeof(double);
    _pointer = std::move(alloc_matrix_in_gpu(size_matrix));
    cudaError_t err =
        cudaMemcpyAsync(_pointer.get(), matrix.data(), size_matrix,
                        cudaMemcpyHostToDevice, stream);
    if (err != 0) {
      throw std::runtime_error("Error copy arrays to device");
    }
  }

  // Unique pointer with custom delete function
  using double_unique_ptr = std::unique_ptr<double, void (*)(double *)>;

 private:
  friend class CudaPipeline;

  // Allocate memory in the GPU for a matrix
  CudaMatrix(long int nrows, long int ncols)
      : _rows{static_cast<int>(nrows)}, _cols{static_cast<int>(ncols)} {
    size_t size_matrix = this->size() * sizeof(double);
    _pointer = std::move(alloc_matrix_in_gpu(size_matrix));
  }

  double_unique_ptr alloc_matrix_in_gpu(size_t size_matrix) const {
    double *dmatrix;
    checkCuda(cudaMalloc(&dmatrix, size_matrix));
    double_unique_ptr dev_ptr(dmatrix,
                              [](double *x) { checkCuda(cudaFree(x)); });
    return dev_ptr;
  }

  double_unique_ptr _pointer{nullptr,
                             [](double *x) { checkCuda(cudaFree(x)); }};
  int _rows;
  int _cols;
};

}  // namespace xtp
}  // namespace votca

#endif
