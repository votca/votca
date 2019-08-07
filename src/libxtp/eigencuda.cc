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
  cublasDestroy(_handle);
  for (auto &p : _allocated) this->fun_free(p.second);
}

template <typename T>
void EigenCuda<T>::fun_alloc(T **x, std::size_t n) const {
  // Allocate memory in the device
  (_pinned) ? cudaMallocHost(x, n) : cudaMalloc(x, n);
}

template <typename T>
void EigenCuda<T>::fun_free(T *x) const {
  // Deallocate memory from the device
  (_pinned) ? cudaFreeHost(x) : cudaFree(x);
};

template <typename T>
void EigenCuda<T>::free_matrix(unsigned id) {
  // Free Array with id from the device
  fun_free(_allocated.at(id));
  _allocated.erase(id);
}

template <typename T>
unsigned EigenCuda<T>::initialize_Matrix(Mat<T> &A, bool copy_to_device) {
  // Copy two matrices to the device

  // size of the Matrices
  std::size_t size_A = A.size() * sizeof(T);

  // Pointer in the device
  T *dA;

  // Allocate either pageable or pinned memory
  fun_alloc(&dA, size_A);

  // Track the allocated variables
  unsigned id = _counter;
  _allocated.emplace(std::make_pair(id, dA));
  _counter += 1;

  // Transfer data to the GPU
  if (copy_to_device) {
    // Pointers at the host
    T *hA = A.data();
    cudaMemcpy(dA, hA, size_A, cudaMemcpyHostToDevice);
  }

  return id;
}

template <typename T>
void EigenCuda<T>::gemm(Shapes sh,
                        std::tuple<unsigned, unsigned, unsigned> ids) {
  // Invoke the gemm subroutine from cublas
  unsigned id_A, id_B, id_C;
  std::tie(id_A, id_B, id_C) = ids;

  // Pointer to the arrays in the device
  T *dA, *dB, *dC;
  dA = _allocated.at(id_A);
  dB = _allocated.at(id_B);
  dC = _allocated.at(id_C);

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

template <typename T>
Mat<T> EigenCuda<T>::dot(Mat<T> &A, Mat<T> &B) {
  // Matrix multiplication

  // Matrix to store the result
  Mat<T> C = Mat<T>::Zero(A.rows(), B.cols());
  std::size_t size_C = C.size() * sizeof(T);

  // Indices of the matrices on the device
  std::tuple<unsigned, unsigned, unsigned> ids = std::make_tuple(
      initialize_Matrix(A), initialize_Matrix(B), initialize_Matrix(C, false));

  // process on GPU
  Shapes sh{A.rows(), A.cols(), B.rows(), B.cols(), C.cols()};
  gemm(sh, ids);

  // send data back to CPU
  T *hC = C.data();
  T *dC = _allocated[std::get<2>(ids)];
  cudaMemcpy(hC, dC, size_C, cudaMemcpyDeviceToHost);

  // Free the result from the device
  unsigned id_C = std::get<2>(ids);
  free_matrix(id_C);

  // create an eigen matrix
  C = Eigen::Map<Mat<T>>(hC, A.rows(), B.cols());

  return C;
}

template <typename T>
std::vector<Mat<T>> EigenCuda<T>::triple_tensor_product(
    Mat<T> &A, Mat<T> &C, std::vector<Mat<T>> &tensor) {
  // Perform the triple matrix multiplication A * matrix * C, for the vector
  // of matrices given by tensor
  std::vector<Mat<T>> rs(tensor.size());

  // Copy Matrix A and B to the device
  unsigned id_A = initialize_Matrix(A);
  unsigned id_C = initialize_Matrix(C);

  // allocate space in device for the temporal matrices
  unsigned size_Y = A.rows() * C.cols() * sizeof(T);
  Mat<T> X = Mat<T>::Zero(A.cols(), C.cols());
  Mat<T> Y = Mat<T>::Zero(A.rows(), C.cols());
  Mat<T> matrix = Mat<T>::Zero(A.cols(), C.rows());

  unsigned id_X = initialize_Matrix(X, false);
  unsigned id_Y = initialize_Matrix(Y, false);
  unsigned id_matrix = initialize_Matrix(matrix, false);

  // Iterate over the tensor Using the previous allocated space in the device
  transform(tensor.begin(), tensor.end(), rs.begin(),
            [this, id_A, id_C, id_X, id_Y, id_matrix, size_Y, &A, &C, &X,
             &Y](Mat<T> &mtx) {
              assert(A.cols() == mtx.rows());
              assert(mtx.cols() == C.rows());

              // Copy matrix to the device
              T *d_matrix = _allocated.at(id_matrix);
              T *h_mtx = mtx.data();

              // move temporal matrix to the preallocated space
              std::size_t size_mtx = mtx.rows() * mtx.cols() * sizeof(T);
              cudaMemcpy(d_matrix, h_mtx, size_mtx, cudaMemcpyHostToDevice);

              // Compute first matrix multiplication
              Shapes sh1{mtx.rows(), mtx.cols(), C.rows(), C.cols(), X.rows()};
              std::tuple<unsigned, unsigned, unsigned> ids =
                  std::make_tuple(id_matrix, id_C, id_X);
              gemm(sh1, ids);

              // compute the second matrix multiplication
              Shapes sh2{A.rows(), A.cols(), X.rows(), X.cols(), Y.rows()};
              ids = std::make_tuple(id_A, id_X, id_Y);
              gemm(sh2, ids);

              // send data back to CPU
              T *hY = Y.data();
              T *dY = this->_allocated[id_Y];
              cudaMemcpy(hY, dY, size_Y, cudaMemcpyDeviceToHost);
              Y = Eigen::Map<Mat<T>>(hY, A.rows(), C.cols());

              return Y;
            });

  return rs;
}

// explicit instantiations
template class EigenCuda<float>;
template class EigenCuda<double>;
}  // namespace xtp
}  // namespace votca
