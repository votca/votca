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

#ifndef VOTCA_XTP_CUDAPIPELINE_H
#define VOTCA_XTP_CUDAPIPELINE_H

// CMake generated file
#include "votca_xtp_config.h"
#ifndef USE_CUDA
#error Cuda not enabled
#endif

#include <type_traits>
// Local VOTCA includes
#include "cudamatrix.h"

/*
 * \brief Perform Tensor-matrix multiplications in a GPU
 *
 * The `CudaPipeline` class handles the allocation and deallocation of arrays on
 * the GPU.
 */

namespace votca {
namespace xtp {

/* \brief The CudaPipeline class offload Eigen operations to an *Nvidia* GPU
 * using the CUDA language. The Cublas handle is the context manager for all the
 * resources needed by Cublas. While a stream is a queue of sequential
 * operations executed in the Nvidia device.
 */
class CudaPipeline {
 public:
  CudaPipeline(int deviceID) : deviceID_{deviceID} {
    checkCuda(cudaSetDevice(deviceID));
    cublasCreate(&handle_);
    cudaStreamCreate(&stream_);
  }

  ~CudaPipeline();
  CudaPipeline() = delete;
  CudaPipeline(const CudaPipeline &) = delete;
  CudaPipeline &operator=(const CudaPipeline &) = delete;

  // C= A*b.asDiagonal()
  template <class M>
  void diag_gemm(const M &A, const CudaMatrix &b, CudaMatrix &C) const;

  // B+=alpha*A;
  void axpy(const CudaMatrix &A, CudaMatrix &B, double alpha = 1.0) const;

  template <class M1, class M2, class M3>
  void gemm(M1 &&A, M2 &&B, M3 &&C, double beta = 0.0) const;

  const cudaStream_t &get_stream() const { return stream_; };

  int getDeviceId() const { return deviceID_; }

 private:
  int deviceID_ = 0;
  // The cublas handles allocates hardware resources on the host and device.
  cublasHandle_t handle_;

  // Asynchronous stream

  cudaStream_t stream_;
};

/*
 * Call the gemm function from cublas, resulting in the multiplication of the
 * two matrices
 */
template <class M1, class M2, class M3>
inline void CudaPipeline::gemm(M1 &&A, M2 &&B, M3 &&C, double beta) const {

  using m1 = std::decay_t<M1>;
  using m2 = std::decay_t<M2>;
  using m3 = std::decay_t<M3>;
  static_assert(!m3::transposed(), "C in gemm cannot be transposed atm");
  // Scalar constanst for calling blas
  double alpha = 1.;
  const double *palpha = &alpha;
  const double *pbeta = &beta;
  cublasOperation_t transA = CUBLAS_OP_N;
  int k = int(A.cols());
  if (m1::transposed()) {
    transA = CUBLAS_OP_T;
    k = int(A.rows());
  }
  cublasOperation_t transB = CUBLAS_OP_N;
  int k2 = int(B.rows());
  if (m2::transposed()) {
    transB = CUBLAS_OP_T;
    k2 = int(B.cols());
  }

  if (k != k2) {
    throw std::runtime_error(
        "Shape mismatch in cuda gemm " + std::to_string(k) + ":" +
        std::to_string(k2) + " A:" + OutputDimension(A) +
        " B:" + OutputDimension(B) + " C:" + OutputDimension(C));
  }

  cublasSetStream(handle_, stream_);
  cublasStatus_t status =
      cublasDgemm(handle_, transA, transB, int(C.rows()), int(C.cols()), k,
                  palpha, A.data(), int(A.ld()), B.data(), int(B.ld()), pbeta,
                  C.data(), int(C.ld()));
  if (status != CUBLAS_STATUS_SUCCESS) {
    throw std::runtime_error("dgemm failed on gpu " +
                             std::to_string(deviceID_) +
                             " with errorcode:" + cudaGetErrorEnum(status));
  }
}

template <class M>
inline void CudaPipeline::diag_gemm(const M &A, const CudaMatrix &b,
                                    CudaMatrix &C) const {

  if (b.cols() != 1 && b.rows() != 1) {
    throw std::runtime_error("B Matrix in Cublas diag_gemm must be a vector");
  }

  cublasSideMode_t mode = CUBLAS_SIDE_RIGHT;
  Index Adim = A.cols();
  if (M::transposed()) {
    mode = CUBLAS_SIDE_LEFT;
    Adim = A.rows();
  }

  if (Adim != b.size()) {
    throw std::runtime_error("Shape mismatch in cuda diag_gemm: A" +
                             OutputDimension(A) + " b" + OutputDimension(b));
  }

  cublasSetStream(handle_, stream_);
  cublasStatus_t status =
      cublasDdgmm(handle_, mode, int(A.rows()), int(A.cols()), A.data(),
                  int(A.ld()), b.data(), 1, C.data(), int(C.ld()));

  if (status != CUBLAS_STATUS_SUCCESS) {
    throw std::runtime_error("diag_gemm failed on gpu " +
                             std::to_string(deviceID_) +
                             " with errorcode:" + cudaGetErrorEnum(status));
  }
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CUDAPIPELINE_H
