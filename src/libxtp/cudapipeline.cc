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

#include <votca/xtp/cudapipeline.h>

namespace votca {
namespace xtp {

CudaPipeline::~CudaPipeline() {

  // destroy handle
  cublasDestroy(_handle);
  // destroy stream
  cudaStreamDestroy(_stream);
}

/*
 * Call the gemm function from cublas, resulting in the multiplication of the
 * two matrices
 */
void CudaPipeline::gemm(const CudaMatrix &A, const CudaMatrix &B,
                        CudaMatrix &C) const {

  // Scalar constanst for calling blas
  double alpha = 1.;
  double beta = 0.;
  const double *palpha = &alpha;
  const double *pbeta = &beta;

  if ((A.cols() != B.rows())) {
    throw std::runtime_error("Shape mismatch in Cublas gemm");
  }
  cublasDgemm(_handle, CUBLAS_OP_N, CUBLAS_OP_N, int(A.rows()), int(B.cols()),
              int(A.cols()), palpha, A.data(), int(A.rows()), B.data(),
              int(B.rows()), pbeta, C.data(), int(C.rows()));
}

}  // namespace xtp
}  // namespace votca
