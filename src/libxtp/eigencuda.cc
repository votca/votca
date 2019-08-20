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
/*
 * Removed all the allocated arrays from the device
 */
template <typename T>
EigenCuda<T>::~EigenCuda() {
  // wait for everything to finish
  cudaDeviceSynchronize();
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
  (_pinned) ? cudaMallocHost(x, n) : cudaMalloc(x, n);
}

/*
 * Deallocate memory from the device
 */
template <typename T>
void EigenCuda<T>::gpu_free(T *x) const {
  (_pinned) ? cudaFreeHost(x) : cudaFree(x);
};

/*
 * Allocate memory in the device for matrix A, then if if `copy_to_device`
 * copy the array to the device. Sometimes it only neccesary to allocate
 * space in the device without copying the array because the initial
 * values may not be important like a temporal matrix.
 */
template <typename T>
T *EigenCuda<T>::initialize_matrix_mem(const Mat<T> &A, bool copy_to_device) {

  // size of the Matrices
  std::size_t size_A = A.size() * sizeof(T);

  // Pointer in the device
  T *dA;

  // Allocate either pageable or pinned memory
  gpu_alloc(&dA, size_A);

  // Transfer data to the GPU
  if (copy_to_device) {
    // Pointers at the host
    const T *hA = A.data();
    cudaError_t err =
        cudaMemcpyAsync(dA, hA, size_A, cudaMemcpyHostToDevice, _stream);
    if (err != 0) {
      throw std::runtime_error("Error copy arrays to device");
    }
  }
  return dA;
}

/*
 * Call the gemm function from cublas, resulting in the multiplication of the
 * two matrices with identifiers id_A and id_B. The result is stored in
 * a Matrix (pointer) with identifier id_C.
 */
template <typename T>
void EigenCuda<T>::gemm(Shapes sh, const T *dA, const T *dB, T *dC) {

  // Scalar constanst for calling blas
  T _alpha = 1.;
  T _beta = 0.;
  const T *_pa = &_alpha;
  const T *_pb = &_beta;

  // call gemm from cublas
  if constexpr (std::is_same<float, T>()) {
    cublasSgemm(_handle, CUBLAS_OP_N, CUBLAS_OP_N, sh.A_rows, sh.B_cols,
                sh.A_cols, _pa, dA, sh.A_rows, dB, sh.B_rows, _pb, dC,
                sh.C_rows);
  } else if (std::is_same<double, T>()) {
    cublasDgemm(_handle, CUBLAS_OP_N, CUBLAS_OP_N, sh.A_rows, sh.B_cols,
                sh.A_cols, _pa, dA, sh.A_rows, dB, sh.B_rows, _pb, dC,
                sh.C_rows);
  }
}

/*
 * Call the gemm?Batched function from cublas, resembling a zip
 * operation of two 3D tensors, applying the gemm function in each pair
 * of matrices. Resulting in a 3D tensor containing the results of the
 * multiplication.
 */
template <typename T>
void EigenCuda<T>::gemmBatched(Shapes sh, const T **dA, const T **dB, T **dC,
                               int batchCount) {

  // Scalar constanst for calling blas
  T _alpha = 1.;
  T _beta = 0.;
  const T *_pa = &_alpha;
  const T *_pb = &_beta;

  // call gemm from cublas
  cublasStatus_t status;
  if constexpr (std::is_same<float, T>()) {
    status = cublasSgemmBatched(_handle, CUBLAS_OP_N, CUBLAS_OP_N, sh.A_rows,
                                sh.B_cols, sh.A_cols, _pa, dA, sh.A_rows, dB,
                                sh.B_rows, _pb, dC, sh.C_rows, batchCount);
  } else if (std::is_same<double, T>()) {
    status = cublasDgemmBatched(_handle, CUBLAS_OP_N, CUBLAS_OP_N, sh.A_rows,
                                sh.B_cols, sh.A_cols, _pa, dA, sh.A_rows, dB,
                                sh.B_rows, _pb, dC, sh.C_rows, batchCount);
  }
  if (status != 0) {
    std::runtime_error("error calling cublas?DgemmBatched");
  }
}

/*
 * \brief Multiply a matrix A by a 3D tensor represented as a vector of
 * matrices. \return a matrix where each column represent the result product.
 * Initially, it allocates memory and copy the matrices A and C together with
 * the tensor to the device. Also, the function allocates the result tensor Y.
 * The method iterates over each submatrix of the tensor computing:
 * Y(i) = tensor(i) * A.
 * Finally, the tensor Y is copy back to the main memory.
 */
template <typename T>
Mat<T> EigenCuda<T>::dot(const Mat<T> &A, const Mat<T> &B) {
  // Matrix to store the result
  Mat<T> C = Mat<T>::Zero(A.rows(), B.cols());
  std::size_t size_C = C.size() * sizeof(T);

  // Indices of the matrices on the device
  T *dA = initialize_matrix_mem(A);
  T *dB = initialize_matrix_mem(B);
  T *dC = initialize_matrix_mem(C, false);

  // process on GPU
  Shapes sh{A.rows(), A.cols(), B.rows(), B.cols(), C.rows()};
  gemm(sh, dA, dB, dC);

  // send data back to CPU
  T *hC = C.data();
  cudaMemcpy(hC, dC, size_C, cudaMemcpyDeviceToHost);

  // create an eigen matrix
  C = Eigen::Map<Mat<T>>(hC, A.rows(), B.cols());

  // Free the result from the device
  gpu_free(dA);
  gpu_free(dB);
  gpu_free(dC);

  return C;
}
template <typename T>
std::vector<Mat<T>> EigenCuda<T>::right_matrix_tensor(
    const Mat<T> &B, const std::vector<Mat<T>> &tensor) {
  // Number of submatrices in the input tensor
  int batchCount = tensor.size();

  // Copy Matrix B to the device
  T *mtxB = initialize_matrix_mem(B);

  // allocate space in device for the temporal matrix
  int rows = tensor[0].rows();  // rows of the submatrices
  int cols = tensor[0].cols();  // cols of the submatrices
  Mat<T> matrix = Mat<T>::Zero(rows, cols);

  // Allocate space and copy to the device the input tensor
  // Notice that hA, hB and hC are arrays IN THE HOST by the pointers
  // are allocated in the DEVICE.
  T *hA[batchCount];
  for (auto i = 0; i < batchCount; i++) {
    hA[i] = initialize_matrix_mem(tensor[i]);
  }

  // represent the matrix B as a tensor where all the submatrices are the same
  T *hB[batchCount];
  for (auto i = 0; i < batchCount; i++) {
    hB[i] = mtxB;
  }

  // Allocate space in the device for the output tensor
  T *hC[batchCount];
  Mat<T> output = Mat<T>::Zero(matrix.rows(), B.cols());
  for (auto i = 0; i < batchCount; i++) {
    hC[i] = initialize_matrix_mem(output, false);
  }

  // Allocate space in the device for the array of pointers
  const T **dA, **dB;
  T **dC;
  size_t size_batch = batchCount * sizeof(T *);
  cudaMalloc(&dA, size_batch);
  cudaMalloc(&dB, size_batch);
  cudaMalloc(&dC, size_batch);

  // Copy the arrays of pointers from host to the device
  cudaMemcpyAsync(dA, hA, size_batch, cudaMemcpyHostToDevice, _stream);
  cudaMemcpyAsync(dB, hB, size_batch, cudaMemcpyHostToDevice, _stream);
  cudaMemcpyAsync(dC, hC, size_batch, cudaMemcpyHostToDevice, _stream);

  // Call tensor matrix multiplication
  Shapes sh{matrix.rows(), matrix.cols(), B.rows(), B.cols(), matrix.rows()};
  gemmBatched(sh, dA, dB, dC, batchCount);

  // Vector containing the results
  std::vector<Mat<T>> rs(batchCount,
                         Mat<T>::Zero(output.rows(), output.cols()));
  std::size_t size_out = output.size() * sizeof(T);

  // Copy Array of pointers on the device to the host
  cudaMemcpyAsync(hC, dC, size_batch, cudaMemcpyDeviceToHost, _stream);

  // Copy each array back to the device
  for (auto i = 0; i < batchCount; i++) {
    T *hout = rs[i].data();
    T *dout = hC[i];
    cudaMemcpyAsync(hout, dout, size_out, cudaMemcpyDeviceToHost, _stream);
    rs[i] = Eigen::Map<Mat<T>>(hout, output.rows(), output.cols());
    ;
  }
  // Deallocate all the memory from the device
  gpu_free(mtxB);
  cudaFree(dA);
  cudaFree(dB);
  cudaFree(dC);
  for (auto i = 0; i < batchCount; i++) {
    gpu_free(hA[i]);
    gpu_free(hC[i]);
  }

  return rs;
}

// explicit instantiations
template class EigenCuda<float>;
template class EigenCuda<double>;

}  // namespace xtp
}  // namespace votca
