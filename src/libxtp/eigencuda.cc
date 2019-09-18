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
 * Allocate memory in the device using either pinned or pageable (default)
 * memory
 */

void EigenCuda::gpu_alloc(double **x, std::size_t n) const {
  (_pinned) ? checkCuda(cudaMallocHost(x, n)) : checkCuda(cudaMalloc(x, n));
}

/*
 * Deallocate memory from the device
 */

void EigenCuda::gpu_free(double *x) const {
  (_pinned) ? checkCuda(cudaFreeHost(x)) : checkCuda(cudaFree(x));
};

/*
 * Check if the available memory is enough to compute the system
 */

void EigenCuda::check_available_memory(size_t requested) const {
  size_t *free, *total;

  // Use Unified memory
  cudaMallocManaged(&free, sizeof(size_t));
  cudaMallocManaged(&total, sizeof(size_t));
  checkCuda(cudaMemGetInfo(free, total));

  std::ostringstream oss;
  oss << "There were requested : " << requested << "bytes int the device\n";
  oss << "Device Free memory (bytes): " << *free
      << "\nDevice total Memory (bytes): " << *total << "\n";

  // Raise an error if there is not enough total or free memory in the device
  if (requested > *free) {
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

double *EigenCuda::initialize_matrix_mem(size_t size_A) const {

  // Pointer in the device
  double *dA;

  // Allocate either pageable or pinned memory
  gpu_alloc(&dA, size_A);

  return dA;
}

/*
 * Allocate memory for the matrix and copy it to the device
 */

double *EigenCuda::initialize_and_copy(const Mat &A) const {

  // allocate memory in the device
  size_t size_A = A.size() * sizeof(double);
  double *dA = initialize_matrix_mem(size_A);

  // Transfer data to the GPU
  const double *hA = A.data();  // Pointers at the host
  cudaError_t err =
      cudaMemcpyAsync(dA, hA, size_A, cudaMemcpyHostToDevice, _stream);
  if (err != 0) {
    throw std::runtime_error("Error copy arrays to device");
  }
  return dA;
}

/*
 * Call the gemm function from cublas, resulting in the multiplication of the
 * two matrices.
 */

void EigenCuda::gemm(Shapes sh, const double *dA, const double *dB,
                     double *dC) const {

  // Scalar constanst for calling blas
  double _alpha = 1.;
  double _beta = 0.;
  const double *_palpha = &_alpha;
  const double *_pbeta = &_beta;

  cublasDgemm(_handle, CUBLAS_OP_N, CUBLAS_OP_N, sh.A_rows, sh.B_cols,
              sh.A_cols, _palpha, dA, sh.A_rows, dB, sh.B_rows, _pbeta, dC,
              sh.C_rows);
}

/*
 * \brief Multiply a matrix B by a 3D tensor represented as a vector of
 * matrices.
 * \return vector of matrices representing the result
 * Initially, it allocates memory and copy the matrix B and each submatrix
 * from tensor to the device. Also, the function allocates the result matrix
 * C into the device. The method iterates
 * over each submatrix of the tensor computing: C = tensor(i) * A.
 */

std::vector<Mat> EigenCuda::right_matrix_tensor(
    const Mat &B, const std::vector<Mat> &tensor) const {
  // Number of submatrices in the input tensor
  int batchCount = tensor.size();

  // First submatrix from the tensor
  Mat matrix = tensor[0];

  // sizes of the matrices to allocated in the device
  size_t size_A = matrix.size() * sizeof(double);
  size_t size_B = B.size() * sizeof(double);
  size_t size_C = matrix.rows() * B.cols() * sizeof(double);

  // Check if there is enough available memory
  check_available_memory(size_A + size_B + size_C);

  // Initialize memory for tensor components
  double *dA = initialize_matrix_mem(size_A);

  // Allocate memory for the final result array C
  double *dC = initialize_matrix_mem(size_C);

  // Copy Matrix B to the device
  double *dB = initialize_and_copy(B);

  // Shapes of the resulting matrices
  Shapes sh{matrix.rows(), matrix.cols(), B.rows(), B.cols(), matrix.rows()};

  // Vector containing the results
  std::vector<Mat> rs(batchCount, Mat::Zero(matrix.rows(), B.cols()));

  // Call tensor matrix multiplication
  for (auto i = 0; i < batchCount; i++) {
    // Copy tensor component to the device
    checkCuda(cudaMemcpyAsync(dA, tensor[i].data(), size_C,
                              cudaMemcpyHostToDevice, _stream));

    // matrix multiplication
    gemm(sh, dA, dB, dC);

    // Copy the result to the host
    double *hout = rs[i].data();
    checkCuda(
        cudaMemcpyAsync(hout, dC, size_C, cudaMemcpyDeviceToHost, _stream));
  }

  // Deallocate all the memory from the device
  checkCuda(cudaFree(dA));
  checkCuda(cudaFree(dB));
  checkCuda(cudaFree(dC));

  return rs;
}

/*
 * \brief performs a matrix_1 * tensor * matrix_2 multiplication
 * \return vector containging the matrix-matrix multiplications
 */
std::vector<Mat> EigenCuda::triple_tensor_product(
    const Mat &A, const Mat &C, const std::vector<Mat> &tensor) {
  // Number of submatrices in the input tensor
  int batchCount = tensor.size();

  // First submatrix from the tensor
  Mat matrix = tensor[0];

  // sizes of the matrices to allocated in the device
  size_t size_A = A.size() * sizeof(double);
  size_t size_B = matrix.size() * sizeof(double);
  size_t size_C = C.size() * sizeof(double);
  std::size_t size_X = A.rows() * matrix.cols() * sizeof(double);
  std::size_t size_Y = A.rows() * C.cols() * sizeof(double);

  // Check if there is enough available memory
  check_available_memory(size_A + size_B + size_C + size_X + size_Y);

  // Copy Matrix B to the device
  double *dA = initialize_and_copy(A);
  double *dC = initialize_and_copy(C);
  double *dB = initialize_matrix_mem(size_B);

  // Intermediate result X
  double *dX = initialize_matrix_mem(size_X);

  // Final result array Y
  double *dY = initialize_matrix_mem(size_Y);

  // Shapes of the matrices
  Shapes sh1{A.rows(), A.cols(), matrix.rows(), matrix.cols(), A.rows()};
  Shapes sh2{A.rows(), matrix.cols(), C.rows(), C.cols(), A.rows()};

  // Vector containing the results
  std::vector<Mat> rs(batchCount, Mat::Zero(A.rows(), C.cols()));

  for (auto i = 0; i < batchCount; i++) {
    // tensor component
    checkCuda(cudaMemcpyAsync(dB, tensor[i].data(), size_B,
                              cudaMemcpyHostToDevice, _stream));

    // Call the first tensor matrix multiplication
    gemm(sh1, dA, dB, dX);

    // Call the second tensor matrix multiplication
    gemm(sh2, dX, dC, dY);

    // Copy the result Array back to the device
    double *hout = rs[i].data();
    checkCuda(
        cudaMemcpyAsync(hout, dY, size_Y, cudaMemcpyDeviceToHost, _stream));
  }

  // Deallocate all the memory from the device
  checkCuda(cudaFree(dA));
  checkCuda(cudaFree(dB));
  checkCuda(cudaFree(dC));
  checkCuda(cudaFree(dX));
  checkCuda(cudaFree(dY));

  return rs;
}
}  // namespace xtp
}  // namespace votca
