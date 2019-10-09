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

void CudaPipeline::throw_if_not_enough_memory_in_gpu(
    size_t requested_memory) const {
  size_t free, total;
  checkCuda(cudaMemGetInfo(&free, &total));

  std::ostringstream oss;
  oss << "There were requested : " << requested_memory
      << "bytes int the device\n";
  oss << "Device Free memory (bytes): " << free
      << "\nDevice total Memory (bytes): " << total << "\n";

  // Raise an error if there is not enough total or free memory in the device
  if (requested_memory > free) {
    oss << "There is not enough memory in the Device!\n";
    throw std::runtime_error(oss.str());
  }
}

/*
 * Call the gemm function from cublas, resulting in the multiplication of the
 * two matrices.
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

  cublasDgemm(_handle, CUBLAS_OP_N, CUBLAS_OP_N, A.rows(), B.cols(), A.cols(),
              palpha, A.pointer(), A.rows(), B.pointer(), B.rows(), pbeta,
              C.pointer(), C.rows());
}

/*
 * \brief Perform a Tensor3D matrix multiplication
 */
void CudaPipeline::right_matrix_tensor_mult(
    std::vector<Eigen::MatrixXd> &tensor, const Eigen::MatrixXd &B) const {
  // First submatrix from the tensor
  const Eigen::MatrixXd &submatrix = tensor[0];

  // sizes of the matrices to allocated in the device
  size_t size_A = submatrix.size() * sizeof(double);
  size_t size_B = B.size() * sizeof(double);
  size_t size_C = submatrix.rows() * B.cols() * sizeof(double);
  throw_if_not_enough_memory_in_gpu(size_A + size_B + size_C);

  // Matrix in the Cuda device
  CudaMatrix matrixA(submatrix.rows(), submatrix.cols());
  CudaMatrix matrixB{B, _stream};
  CudaMatrix matrixC(submatrix.rows(), B.cols());

  // Call tensor matrix multiplication
  for (auto i = 0; i < static_cast<int>(tensor.size()); i++) {
    matrixA.copy_to_gpu(tensor[i]);
    gemm(matrixA, matrixB, matrixC);
    // Copy the result to the host
    tensor[i] = matrixC;
  }
}

/*
 * \brief performs a CudaMatrix * EigenMatrix * CudaMatrix multiplication
 */
Eigen::MatrixXd CudaPipeline::triple_matrix_mult(const CudaMatrix &A,
                                                 const Eigen::MatrixXd &matrix,
                                                 const CudaMatrix &C) const {
  // sizes of the matrices to allocated in the device
  size_t size_B = matrix.size() * sizeof(double);
  std::size_t size_W = A.rows() * matrix.cols() * sizeof(double);
  std::size_t size_Z = A.rows() * C.cols() * sizeof(double);
  throw_if_not_enough_memory_in_gpu(size_B + size_W + size_Z);

  // Intermediate Matrices
  CudaMatrix B{matrix, _stream};
  CudaMatrix W{A.rows(), matrix.cols()};
  CudaMatrix Z{A.rows(), C.cols()};
  gemm(A, B, W);
  gemm(W, C, Z);

  // Copy the result Array back to the device
  Eigen::MatrixXd result = Z;

  return result;
}

}  // namespace xtp
}  // namespace votca
