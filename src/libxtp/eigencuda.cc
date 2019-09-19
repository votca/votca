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
 * Deallocate memory from the device
 */
void free_mem_in_gpu(double *x) { checkCuda(cudaFree(x)); };

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
EigenCuda::uniq_double EigenCuda::alloc_matrix_in_gpu(
    size_t size_matrix) const {

  // Pointer in the device
  double *dmatrix;
  checkCuda(cudaMalloc(&dmatrix, size_matrix));
  uniq_double dev_ptr(dmatrix, free_mem_in_gpu);
  return dev_ptr;
}

EigenCuda::uniq_double EigenCuda::copy_matrix_to_gpu(const Mat &matrix) const {
  // allocate memory in the device
  size_t size_matrix = matrix.size() * sizeof(double);
  uniq_double dev_ptr = alloc_matrix_in_gpu(size_matrix);

  // Transfer data to the GPU
  const double *hmatrix = matrix.data();  // Pointers at the host
  cudaError_t err = cudaMemcpyAsync(dev_ptr.get(), hmatrix, size_matrix,
                                    cudaMemcpyHostToDevice, _stream);
  if (err != 0) {
    throw std::runtime_error("Error copy arrays to device");
  }
  return dev_ptr;
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
  check_available_memory_in_gpu(size_A + size_B + size_C);

  uniq_double dA = alloc_matrix_in_gpu(size_A);
  uniq_double dC = alloc_matrix_in_gpu(size_C);
  uniq_double dB = copy_matrix_to_gpu(B);

  ShapesOfMatrices sh{matrix.rows(), matrix.cols(), B.rows(), B.cols(),
                      matrix.rows()};

  std::vector<Mat> result(batchCount, Mat::Zero(matrix.rows(), B.cols()));

  // Call tensor matrix multiplication
  for (auto i = 0; i < batchCount; i++) {
    // Copy tensor component to the device
    checkCuda(cudaMemcpyAsync(dA.get(), tensor[i].data(), size_C,
                              cudaMemcpyHostToDevice, _stream));

    // matrix multiplication
    gemm(sh, dA.get(), dB.get(), dC.get());

    // Copy the result to the host
    double *hout = result[i].data();
    checkCuda(cudaMemcpyAsync(hout, dC.get(), size_C, cudaMemcpyDeviceToHost,
                              _stream));
  }

  return result;
}

}  // namespace xtp
}  // namespace votca
