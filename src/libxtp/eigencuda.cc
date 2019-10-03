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

#include <votca/xtp/eigencuda.h>

namespace votca {
namespace xtp {

EigenCuda::~EigenCuda() {

  // destroy handle
  cublasDestroy(_handle);
  // destroy stream
  cudaStreamDestroy(_stream);
}

/*
 * Check if the available memory is enough to compute the system
 */
void EigenCuda::check_available_memory_in_gpu(size_t requested_memory) const {
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
 * Allocate memory in the device for matrix A.
 */
CudaMatrix::double_unique_ptr CudaMatrix::alloc_matrix_in_gpu(
    size_t size_matrix) const {

  // Pointer in the device
  double *dmatrix;
  checkCuda(cudaMalloc(&dmatrix, size_matrix));
  double_unique_ptr dev_ptr(dmatrix, [](double *x) { checkCuda(cudaFree(x)); });
  return dev_ptr;
}

void CudaMatrix::copy_matrix_to_gpu(const Eigen::MatrixXd &matrix,
                                    const cudaStream_t &stream) const {
  // Transfer data to the GPU
  size_t size_matrix = matrix.size() * sizeof(double);
  const double *hmatrix = matrix.data();  // Pointers at the host
  cudaError_t err = cudaMemcpyAsync(_pointer.get(), hmatrix, size_matrix,
                                    cudaMemcpyHostToDevice, stream);
  if (err != 0) {
    throw std::runtime_error("Error copy arrays to device");
  }
}

/*
 * Call the gemm function from cublas, resulting in the multiplication of the
 * two matrices.
 */
void EigenCuda::gemm(const CudaMatrix &A, const CudaMatrix &B,
                     CudaMatrix &C) const {

  // Scalar constanst for calling blas
  double alpha = 1.;
  double beta = 0.;
  const double *palpha = &alpha;
  const double *pbeta = &beta;

  cublasDgemm(_handle, CUBLAS_OP_N, CUBLAS_OP_N, A.rows(), B.cols(), A.cols(),
              palpha, A.pointer(), A.rows(), B.pointer(), B.rows(), pbeta,
              C.pointer(), C.rows());
}

/*
 * \brief Perform a Tensor3D matrix multiplication
 */
void EigenCuda::right_matrix_tensor_mult(std::vector<Eigen::MatrixXd> &tensor,
                                         const Eigen::MatrixXd &B) const {
  // First submatrix from the tensor
  const Eigen::MatrixXd &submatrix = tensor[0];

  // sizes of the matrices to allocated in the device
  size_t size_A = submatrix.size() * sizeof(double);
  size_t size_B = B.size() * sizeof(double);
  size_t size_C = submatrix.rows() * B.cols() * sizeof(double);
  check_available_memory_in_gpu(size_A + size_B + size_C);

  // Matrix in the Cuda device

  CudaMatrix matrixA{submatrix.rows(), submatrix.cols()};
  CudaMatrix matrixB{B.rows(), B.cols()};
  CudaMatrix matrixC{submatrix.rows(), B.cols()};
  matrixB.copy_matrix_to_gpu(B, _stream);

  // Call tensor matrix multiplication
  for (auto i = 0; i < static_cast<int>(tensor.size()); i++) {
    // Copy tensor component to the device
    checkCuda(cudaMemcpyAsync(matrixA.pointer(), tensor[i].data(), size_C,
                              cudaMemcpyHostToDevice, _stream));

    // matrix multiplication
    gemm(matrixA, matrixB, matrixC);

    // Copy the result to the host
    double *hout = tensor[i].data();
    checkCuda(cudaMemcpyAsync(hout, matrixC.pointer(), size_C,
                              cudaMemcpyDeviceToHost, _stream));
  }
}

/*
 * \brief performs a matrix_1 * matrix2 * matrix_2 multiplication
 */
Eigen::MatrixXd EigenCuda::triple_matrix_mult(const CudaMatrix &A,
                                              const Eigen::MatrixXd &matrix,
                                              const CudaMatrix &C) const {

  // sizes of the matrices to allocated in the device
  size_t size_B = matrix.size() * sizeof(double);
  std::size_t size_W = A.rows() * matrix.cols() * sizeof(double);
  std::size_t size_Z = A.rows() * C.cols() * sizeof(double);

  // Check if there is enough available memory
  check_available_memory_in_gpu(size_B + size_W + size_Z);

  // Intermediate Matrices
  CudaMatrix B{matrix.rows(), matrix.cols()};
  CudaMatrix W{A.rows(), matrix.cols()};
  CudaMatrix Z{A.rows(), C.cols()};
  B.copy_matrix_to_gpu(matrix, _stream);

  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(A.rows(), C.cols());

  // Call the first tensor matrix multiplication
  gemm(A, B, W);

  // Call the second tensor matrix multiplication
  gemm(W, C, Z);

  // Copy the result Array back to the device
  checkCuda(cudaMemcpyAsync(result.data(), Z.pointer(), size_Z,
                            cudaMemcpyDeviceToHost, _stream));

  return result;
}

}  // namespace xtp
}  // namespace votca
