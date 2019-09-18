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
 * Allocate memory in the gpu using either pinned or pageable (default)
 * memory
 */
void EigenCuda::alloc_mem_in_gpu(double **x, std::size_t n) const {
  (_pinned) ? checkCuda(cudaMallocHost(x, n)) : checkCuda(cudaMalloc(x, n));
}

/*
 * Deallocate memory from the device
 */
void EigenCuda::free_mem_in_gpu(double *x) const {
  (_pinned) ? checkCuda(cudaFreeHost(x)) : checkCuda(cudaFree(x));
};

/*
 * Check if the available memory is enough to compute the system
 */
void EigenCuda::check_available_memory_in_gpu(size_t requested_memory) const {
  // Use Unified memory
  size_t *free, *total;
  cudaMallocManaged(&free, sizeof(size_t));
  cudaMallocManaged(&total, sizeof(size_t));
  checkCuda(cudaMemGetInfo(free, total));

  std::ostringstream oss;
  oss << "There were requested : " << requested_memory
      << "bytes int the device\n";
  oss << "Device Free memory (bytes): " << *free
      << "\nDevice total Memory (bytes): " << *total << "\n";

  // Raise an error if there is not enough total or free memory in the device
  if (requested_memory > *free) {
    oss << "There is not enough memory in the Device!\n";
    throw std::runtime_error(oss.str());
  }

  // Release memory
  cudaFree(free);
  cudaFree(total);
}

/*
 * Allocate memory in the device for matrix A.
 */
double *EigenCuda::alloc_matrix_in_gpu(size_t size_matrix) const {

  // Pointer in the device
  double *dmatrix;

  // Allocate either pageable or pinned memory
  alloc_mem_in_gpu(&dmatrix, size_matrix);

  return dmatrix;
}

double *EigenCuda::copy_matrix_to_gpu(const Mat &matrix) const {
  // allocate memory in the device
  size_t size_matrix = matrix.size() * sizeof(double);
  double *dmatrix = alloc_matrix_in_gpu(size_matrix);

  // Transfer data to the GPU
  const double *hmatrix = matrix.data();  // Pointers at the host
  cudaError_t err = cudaMemcpyAsync(dmatrix, hmatrix, size_matrix,
                                    cudaMemcpyHostToDevice, _stream);
  if (err != 0) {
    throw std::runtime_error("Error copy arrays to device");
  }
  return dmatrix;
}

/*
 * Call the gemm function from cublas, resulting in the multiplication of the
 * two matrices.
 */
void EigenCuda::gemm(ShapesOfMatrices sh, const double *dA, const double *dB,
                     double *dC) const {

  // Scalar constanst for calling blas
  double alpha = 1.;
  double beta = 0.;
  const double *palpha = &alpha;
  const double *pbeta = &beta;

  cublasDgemm(_handle, CUBLAS_OP_N, CUBLAS_OP_N, sh.A_rows, sh.B_cols,
              sh.A_cols, palpha, dA, sh.A_rows, dB, sh.B_rows, pbeta, dC,
              sh.C_rows);
}

/*
 * \brief Perform a Tensor3D matrix multiplication
 */
std::vector<Mat> EigenCuda::right_matrix_tensor_mult(
    const std::vector<Mat> &tensor, const Mat &B) const {
  int batchCount = tensor.size();

  // First submatrix from the tensor
  const Mat &matrix = tensor[0];

  // sizes of the matrices to allocated in the device
  size_t size_A = matrix.size() * sizeof(double);
  size_t size_B = B.size() * sizeof(double);
  size_t size_C = matrix.rows() * B.cols() * sizeof(double);

  // Check if there is enough available memory
  check_available_memory_in_gpu(size_A + size_B + size_C);

  double *dA = alloc_matrix_in_gpu(size_A);
  double *dC = alloc_matrix_in_gpu(size_C);
  double *dB = copy_matrix_to_gpu(B);

  ShapesOfMatrices sh{matrix.rows(), matrix.cols(), B.rows(), B.cols(),
                      matrix.rows()};

  std::vector<Mat> result(batchCount, Mat::Zero(matrix.rows(), B.cols()));

  // Call tensor matrix multiplication
  for (auto i = 0; i < batchCount; i++) {
    // Copy tensor component to the device
    checkCuda(cudaMemcpyAsync(dA, tensor[i].data(), size_C,
                              cudaMemcpyHostToDevice, _stream));

    // matrix multiplication
    gemm(sh, dA, dB, dC);

    // Copy the result to the host
    double *hout = result[i].data();
    checkCuda(
        cudaMemcpyAsync(hout, dC, size_C, cudaMemcpyDeviceToHost, _stream));
  }

  // Deallocate all the memory from the device
  free_mem_in_gpu(dA);
  free_mem_in_gpu(dB);
  free_mem_in_gpu(dC);

  return result;
}

/*
 * \brief performs a matrix_1 * tensor * matrix_2 multiplication
 */
std::vector<Mat> EigenCuda::matrix_tensor_matrix_mult(
    const Mat &A, const std::vector<Mat> &tensor, const Mat &C) {
  int batchCount = tensor.size();

  // First submatrix from the tensor
  const Mat &matrix = tensor[0];

  // sizes of the matrices to allocated in the device
  size_t size_A = A.size() * sizeof(double);
  size_t size_B = matrix.size() * sizeof(double);
  size_t size_C = C.size() * sizeof(double);
  std::size_t size_X = A.rows() * matrix.cols() * sizeof(double);
  std::size_t size_Y = A.rows() * C.cols() * sizeof(double);

  // Check if there is enough available memory
  check_available_memory_in_gpu(size_A + size_B + size_C + size_X + size_Y);

  double *dA = copy_matrix_to_gpu(A);
  double *dC = copy_matrix_to_gpu(C);
  double *dB = alloc_matrix_in_gpu(size_B);

  // Intermediate result X
  double *dX = alloc_matrix_in_gpu(size_X);

  // Final result array Y
  double *dY = alloc_matrix_in_gpu(size_Y);

  ShapesOfMatrices sh1{A.rows(), A.cols(), matrix.rows(), matrix.cols(),
                       A.rows()};
  ShapesOfMatrices sh2{A.rows(), matrix.cols(), C.rows(), C.cols(), A.rows()};

  std::vector<Mat> result(batchCount, Mat::Zero(A.rows(), C.cols()));
  for (auto i = 0; i < batchCount; i++) {
    // tensor component
    checkCuda(cudaMemcpyAsync(dB, tensor[i].data(), size_B,
                              cudaMemcpyHostToDevice, _stream));

    // Call the first tensor matrix multiplication
    gemm(sh1, dA, dB, dX);

    // Call the second tensor matrix multiplication
    gemm(sh2, dX, dC, dY);

    // Copy the result Array back to the device
    double *hout = result[i].data();
    checkCuda(
        cudaMemcpyAsync(hout, dY, size_Y, cudaMemcpyDeviceToHost, _stream));
  }

  // Deallocate all the memory from the device
  free_mem_in_gpu(dA);
  free_mem_in_gpu(dB);
  free_mem_in_gpu(dC);
  free_mem_in_gpu(dX);
  free_mem_in_gpu(dY);

  return result;
}
}  // namespace xtp
}  // namespace votca
