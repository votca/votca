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

#ifndef __XTP_CUDA_PIPELINE__H
#define __XTP_CUDA_PIPELINE__H

#include <cublas_v2.h>
#include <curand.h>
#include <vector>
#include <votca/xtp/cudamatrix.h>
#include <votca/xtp/eigen.h>

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
  CudaPipeline() {
    cublasCreate(&_handle);
    cudaStreamCreate(&_stream);
  }
  ~CudaPipeline();

  CudaPipeline(const CudaPipeline &) = delete;
  CudaPipeline &operator=(const CudaPipeline &) = delete;

  // Invoke the ?gemm function of cublas
  void gemm(const CudaMatrix &A, const CudaMatrix &B, CudaMatrix &C) const;

  const cudaStream_t &get_stream() const { return _stream; };

 private:
  // The cublas handles allocates hardware resources on the host and device.
  cublasHandle_t _handle;

  // Asynchronous stream
  cudaStream_t _stream;
};

}  // namespace xtp
}  // namespace votca

#endif
