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

template <typename T>
EigenCuda<T>::~EigenCuda() {

  // destroy handle
  cublasDestroy(_handle);
  // destroy stream
  cudaStreamDestroy(_stream);
}

/*
 * Allocate memory in the device using either pinned or pageable (default)
 * memory
 */
template <typename T>
void EigenCuda<T>::gpu_alloc(T **x, std::size_t n) const {
  (_pinned) ? checkCuda(cudaMallocHost(x, n)) : checkCuda(cudaMalloc(x, n));
}

/*
 * Deallocate memory from the device
 */
template <typename T>
void EigenCuda<T>::gpu_free(T *x) const {
  (_pinned) ? checkCuda(cudaFreeHost(x)) : checkCuda(cudaFree(x));
};

/*
 * Allocate memory in the device for matrix A.
 */
template <typename T>
T *EigenCuda<T>::initialize_matrix_mem(size_t size_A) const {

  // Pointer in the device
  T *dA;

  // Allocate either pageable or pinned memory
  gpu_alloc(&dA, size_A);

  return dA;
}

/*
 * Allocate memory for the matrix and copy it to the device
 */
template <typename T>
T *EigenCuda<T>::initialize_and_copy(const Mat<T> &A) const {

  // allocate memory in the device
  size_t size_A = A.size() * sizeof(T);
  T *dA = initialize_matrix_mem(size_A);

  // Transfer data to the GPU 
  const T *hA = A.data();   // Pointers at the host
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
template <typename T>
void EigenCuda<T>::gemm(Shapes sh, const T *dA, const T *dB, T *dC) const {

  // Scalar constanst for calling blas
  T _alpha = 1.;
  T _beta = 0.;
  const T *_palpha = &_alpha;
  const T *_pbeta = &_beta;
  
  // call gemm from cublas
  if constexpr (std::is_same<float, T>()) {
    cublasSgemm(_handle, CUBLAS_OP_N, CUBLAS_OP_N, sh.A_rows, sh.B_cols,
                sh.A_cols, _palpha, dA, sh.A_rows, dB, sh.B_rows, _pbeta, dC,
                sh.C_rows);
  } else if (std::is_same<double, T>()) {
    cublasDgemm(_handle, CUBLAS_OP_N, CUBLAS_OP_N, sh.A_rows, sh.B_cols,
                sh.A_cols, _palpha, dA, sh.A_rows, dB, sh.B_rows, _pbeta, dC,
                sh.C_rows);
  }
}

/*
 * \brief Matrix-Matrix multiplication in GPU
 */
template <typename T>
Mat<T> EigenCuda<T>::dot(const Mat<T> &A, const Mat<T> &B) const {
  // Matrix to store the result
  Mat<T> C = Mat<T>::Zero(A.rows(), B.cols());
  std::size_t size_C = C.size() * sizeof(T);

  // Indices of the matrices on the device
  T *dA = initialize_and_copy(A);
  T *dB = initialize_and_copy(B);
  T *dC = initialize_matrix_mem(size_C);

  // process on GPU
  Shapes sh{A.rows(), A.cols(), B.rows(), B.cols(), C.rows()};
  gemm(sh, dA, dB, dC);

  // send data back to CPU
  T *hC = C.data();
  cudaMemcpy(hC, dC, size_C, cudaMemcpyDeviceToHost);

  // Free the result from the device
  gpu_free(dA);
  gpu_free(dB);
  gpu_free(dC);

  return C;
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
template <typename T>
std::vector<Mat<T>> EigenCuda<T>::right_matrix_tensor(
    const Mat<T> &B, const std::vector<Mat<T>> &tensor) const {
  // Number of submatrices in the input tensor
  int batchCount = tensor.size();

  // First submatrix from the tensor
  Mat<T> matrix = tensor[0];

  // Initialize memory for tensor components
  size_t size_A = matrix.size() * sizeof(T);
  T *dA = initialize_matrix_mem(size_A);
  
  // Copy Matrix B to the device
  T *dB = initialize_and_copy(B);

  // Allocate memory for the final result array C
  size_t size_C = matrix.rows() * B.cols() * sizeof(T);
  T *dC = initialize_matrix_mem(size_C);

  // Shapes of the resulting matrices
  Shapes sh{matrix.rows(), matrix.cols(), B.rows(), B.cols(), matrix.rows()};

  // Vector containing the results
  std::vector<Mat<T>> rs(batchCount, Mat<T>::Zero(matrix.rows(), B.cols()));
  
  // Call tensor matrix multiplication
  for (auto i=0; i < batchCount; i++) {
    // Copy tensor component to the device
    checkCuda(cudaMemcpyAsync(dA, tensor[i].data(), size_C, cudaMemcpyHostToDevice, _stream));
    
    // matrix multiplication
    gemm(sh, dA, dB, dC);

    // Copy the result to the host
    T *hout = rs[i].data();
    checkCuda(cudaMemcpyAsync(hout, dC, size_C, cudaMemcpyDeviceToHost, _stream));
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
template <typename T>
std::vector<Mat<T>> EigenCuda<T>::triple_tensor_product(
    const Mat<T> &A, const Mat<T> &C, const std::vector<Mat<T>> &tensor) {
  // Number of submatrices in the input tensor
  int batchCount = tensor.size();

  // Copy Matrix B to the device
  T *dA = initialize_and_copy(A);
  T *dC = initialize_and_copy(C);

  // First submatrix from the tensor
  Mat<T> matrix = tensor[0];
  size_t size_B = matrix.size() * sizeof(T);
  T *dB = initialize_matrix_mem(size_B);

  // Intermediate result X
  std::size_t size_X = A.rows() * matrix.cols() * sizeof(T);
  T *dX = initialize_matrix_mem(size_X) ;
  
  // Final result array Y
  std::size_t size_Y = A.rows() * C.cols() * sizeof(T);
  T *dY = initialize_matrix_mem(size_Y);

  // Shapes of the matrices
  Shapes sh1{A.rows(), A.cols(), matrix.rows(), matrix.cols(), A.rows()};
  Shapes sh2{A.rows(), matrix.cols(), C.rows(), C.cols(), A.rows()};  

  // Vector containing the results
  std::vector<Mat<T>> rs(batchCount, Mat<T>::Zero(A.rows(), C.cols()));

  for (auto i=0; i < batchCount; i++) {
    // tensor component
    checkCuda(cudaMemcpyAsync(dB, tensor[i].data(), size_B, cudaMemcpyHostToDevice, _stream));

    // Call the first tensor matrix multiplication
    gemm(sh1, dA, dB, dX);

    // Call the second tensor matrix multiplication
    gemm(sh2, dX, dC, dY);

    // Copy the result Array back to the device
    T *hout = rs[i].data();
    checkCuda(cudaMemcpyAsync(hout, dY, size_Y, cudaMemcpyDeviceToHost, _stream));
  }

  // Deallocate all the memory from the device
  checkCuda(cudaFree(dA));
  checkCuda(cudaFree(dB));
  checkCuda(cudaFree(dC));
  checkCuda(cudaFree(dX));
  checkCuda(cudaFree(dY));

  return rs;
}

// explicit instantiations
template class EigenCuda<float>;
template class EigenCuda<double>;

}  // namespace xtp
}  // namespace votca
