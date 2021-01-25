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

namespace votca {
namespace xtp {

CudaPipeline::~CudaPipeline() {

  // destroy handle
  cublasDestroy(_handle);
  // destroy stream
  cudaStreamDestroy(_stream);
}

std::string CudaPipeline::cudaGetErrorEnum(cublasStatus_t error) {
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

void CudaPipeline::diag_gemm(const CudaMatrix &A, const CudaMatrix &b,
                             CudaMatrix &C) const {

  if (b.cols() != 1) {
    throw std::runtime_error(
        "B Matrix in Cublas diag_gemm must have one column");
  }

  if (A.rows() != b.rows()) {
    throw std::runtime_error("Shape mismatch in cuda diag_gemm");
  }

  cublasSetStream(_handle, _stream);
  cublasStatus_t status = cublasDdgmm(_handle, CUBLAS_SIDE_LEFT, int(A.rows()),
                                      int(A.cols()), A.data(), int(A.rows()),
                                      b.data(), 1, C.data(), int(C.rows()));

  if (status != CUBLAS_STATUS_SUCCESS) {
    throw std::runtime_error("diag_gemm failed on gpu " +
                             std::to_string(_deviceID) +
                             " with errorcode:" + cudaGetErrorEnum(status));
  }
}

}  // namespace xtp
}  // namespace votca
