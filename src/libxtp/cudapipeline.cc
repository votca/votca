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
#include "votca/xtp/cudapipeline.h"
#include <stdexcept>
#include <string>

namespace votca {
namespace xtp {

CudaPipeline::~CudaPipeline() {

  // destroy handle
  cublasDestroy(_handle);
  // destroy stream
  cudaStreamDestroy(_stream);
}

void CudaPipeline::axpy(const CudaMatrix &A, CudaMatrix &B,
                        double alpha) const {

  if (A.rows() != B.rows() || A.cols() != B.cols()) {
    throw std::runtime_error("Shape mismatch in cuda axpy");
  }

  cublasSetStream(_handle, _stream);
  cublasStatus_t status =
      cublasDaxpy(_handle, int(A.size()), &alpha, A.data(), 1, B.data(), 1);

  if (status != CUBLAS_STATUS_SUCCESS) {
    throw std::runtime_error("axpy failed on gpu " + std::to_string(_deviceID) +
                             " with errorcode:" + cudaGetErrorEnum(status));
  }
}



}  // namespace xtp
}  // namespace votca
