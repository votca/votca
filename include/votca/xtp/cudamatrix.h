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

#ifndef VOTCA_XTP_CUDAMATRIX_H
#define VOTCA_XTP_CUDAMATRIX_H

// CMake generated file
#include "votca_xtp_config.h"
#include <stdexcept>
#include <string>
#ifndef USE_CUDA
#error Cuda not enabled
#endif

// Standard includes
#include <cassert>
#include <memory>
// Third party includes
#include <cublas_v2.h>
#include <curand.h>

// Local VOTCA includes
#include "eigen.h"

/*
 * \brief Matrix Representation inside an Nvidia GPU
 */
namespace votca {
namespace xtp {

void checkCuda(cudaError_t result);
void checkCublas(cublasStatus_t result);
std::string cudaGetErrorEnum(cublasStatus_t error);
Index count_available_gpus();

template <class M>
std::string OutputDimension(const M &mat) {
  std::string transpose=M::transposed()? "T":"";

  return std::string(transpose+"(" + std::to_string(mat.rows()) + "x" +
                     std::to_string(mat.cols()) + ")");
}

template <class M>
class CudaMatrixBlock {
 public:
  CudaMatrixBlock(const M &mat, Index rowoffset, Index coloffset, Index rows,
                  Index cols)
      : mat_(mat), rows_(rows), cols_(cols) {

    assert((rowoffset + rows) <= mat.rows() &&
           "block has to fit in matrix, rows exceeded");
    assert((coloffset + cols) <= mat.cols() &&
           "block has to fit in matrix, cols exceeded");
    start_ = coloffset * ld() + rowoffset;
  }
  Index size() const { return rows() * cols(); }
  Index rows() const { return rows_; }
  Index cols() const { return cols_; }
  Index ld() const { return mat_.ld(); }
  double *data() const { return mat_.data() + start_; }

  static constexpr bool transposed() { return M::transposed(); }

 private:
  const M &mat_;
  Index rows_;
  Index cols_;
  Index start_;
};

template <class M>
class CudaMatrixTranspose {
 public:
  CudaMatrixTranspose(const M &mat) : mat_(mat) { ; }
  Index size() const { return mat_.size(); }
  Index rows() const { return mat_.rows(); }
  Index cols() const { return mat_.cols(); }
  Index ld() const { return mat_.ld(); }
  double *data() const { return mat_.data(); }

  static constexpr bool transposed() { return !M::transposed(); }

 private:
  const M &mat_;
};

class CudaMatrix {
 public:
  Index size() const { return _ld * _cols; };
  Index rows() const { return _ld; };
  Index ld() const { return _ld; }
  Index cols() const { return _cols; };
  double *data() const { return _data.get(); };

  void reshape(Index rows, Index cols) {
    assert(rows * cols == size() &&
           "reshape cannot change the shape of the matrix");
    _cols = cols;
    _ld = rows;
  }

  CudaMatrixTranspose<CudaMatrix> transpose() const {
    return CudaMatrixTranspose<CudaMatrix>(*this);
  }

  CudaMatrixBlock<CudaMatrix> block(Index rowoffset, Index coloffset,
                                    Index rows, Index cols) const {
    return CudaMatrixBlock<CudaMatrix>(*this, rowoffset, coloffset, rows, cols);
  }

  static constexpr bool transposed() { return false; }

  template <class T>
  void copy_to_gpu(const T &m) {
    if(m.rows()!=_ld || m.cols()!=_cols){
      throw std::runtime_error("Shape mismatch of cpu ("+std::to_string(m.rows())+"x"+std::to_string(m.cols())+") and gpu matrix"+OutputDimension(*this));
    }
    checkCublas(cublasSetMatrixAsync(
        int(m.rows()), int(m.cols()), sizeof(double), m.data(),
        int(m.colStride()), this->data(), int(this->rows()), _stream));
  }

  template <class T>
  CudaMatrix(const T &matrix, const cudaStream_t &stream)
      : _ld{static_cast<Index>(matrix.rows())},
        _cols{static_cast<Index>(matrix.cols())} {
    _data = alloc_matrix_in_gpu(size_matrix());
    _stream = stream;
    copy_to_gpu(matrix);
  }

  // Allocate memory in the GPU for a matrix
  CudaMatrix(Index nrows, Index ncols, const cudaStream_t &stream);

  void setZero();

  // Convert A Cudamatrix to an EigenMatrix
  operator Eigen::MatrixXd() const;

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
  Index _ld;
  Index _cols;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CUDAMATRIX_H
