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
  CudaPipeline(int deviceID) : _deviceID{deviceID} {
    checkCuda(cudaSetDevice(deviceID));
    cublasCreate(&_handle);
    cudaStreamCreate(&_stream);
  }

  ~CudaPipeline();
  CudaPipeline() = delete;
  CudaPipeline(const CudaPipeline &) = delete;
  CudaPipeline &operator=(const CudaPipeline &) = delete;

  // Invoke the ?gemm function of cublas
  template<class M1,class M2>
  void gemm(const M1 &A, const M2 &B, CudaMatrix &C,
            double beta = 0.0) const;

  // Invoke the multiplication with a diagonal matrix of cublas, diagonal matrix
  // B must have 1 column
  void diag_gemm(const CudaMatrix &A, const CudaMatrix &b, CudaMatrix &C) const;

  const cudaStream_t &get_stream() const { return _stream; };

  int getDeviceId() const { return _deviceID; }

 private:
  int _deviceID = 0;
  // The cublas handles allocates hardware resources on the host and device.
  cublasHandle_t _handle;

  // Asynchronous stream
  cudaStream_t _stream;

  static std::string cudaGetErrorEnum(cublasStatus_t error);
};

/*
 * Call the gemm function from cublas, resulting in the multiplication of the
 * two matrices
 */
template<class M1,class M2>
inline void CudaPipeline::gemm(const M1 &A, const M2 &B, CudaMatrix &C, double beta) const {

  // Scalar constanst for calling blas
  double alpha = 1.;
  const double *palpha = &alpha;
  const double *pbeta = &beta;
  cublasOperation_t transA = CUBLAS_OP_N;
  int k = int(A.cols());
  if (A.transposed()) {
    transA = CUBLAS_OP_T;
    k = int(A.rows());
  }
  cublasOperation_t transB = CUBLAS_OP_N;
  int k2 = int(B.rows());
  if (B.transposed()) {
    transB = CUBLAS_OP_T;
    k2 = int(B.cols());
  }

  if (k != k2) {
    throw std::runtime_error("Shape mismatch in cuda gemm");
  }

  cublasSetStream(_handle, _stream);
  cublasStatus_t status =
      cublasDgemm(_handle, transA, transB, int(C.rows()), int(C.cols()), k,
                  palpha, A.data(), int(A.rows()), B.data(), int(B.rows()),
                  pbeta, C.data(), int(C.rows()));
  if (status != CUBLAS_STATUS_SUCCESS) {
    throw std::runtime_error("dgemm failed on gpu " +
                             std::to_string(_deviceID) +
                             " with errorcode:" + cudaGetErrorEnum(status));
  }
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CUDAPIPELINE_H
