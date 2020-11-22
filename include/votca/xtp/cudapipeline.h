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
    cublasCreate(&_handle);
    cudaStreamCreate(&_stream);
    checkCuda(cudaSetDevice(deviceID));
  }

  ~CudaPipeline();
  CudaPipeline() = delete;
  CudaPipeline(const CudaPipeline &) = delete;
  CudaPipeline &operator=(const CudaPipeline &) = delete;

CudaPipeline( CudaPipeline &&) = default;
  // Invoke the ?gemm function of cublas
  void gemm(const CudaMatrix &A, const CudaMatrix &B, CudaMatrix &C,
            bool transpose_A = false, bool transpose_B = false,
            double beta = 0.0) const;

  // Invoke the multiplication with a diagonal matrix of cublas, diagonal matrix
  // B must have 1 column
  void diag_gemm(const CudaMatrix &A, const CudaMatrix &b, CudaMatrix &C) const;

  const cudaStream_t &get_stream() const { return _stream; };

  int getDeviceId()const{return _deviceID;}

 private:
  int _deviceID = 0;
  // The cublas handles allocates hardware resources on the host and device.
  cublasHandle_t _handle;

  // Asynchronous stream
  cudaStream_t _stream;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CUDAPIPELINE_H
