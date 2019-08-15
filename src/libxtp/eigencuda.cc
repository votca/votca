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
 * Stack a vector of matrices as a single matrix, where each row corresponds
 * to a matrix.
 */
template <typename T> Mat<T> stack(const std::vector<Mat<T>> &tensor) {

  int rows = tensor.size();
  int cols = tensor[0].size(); // size of each matrix

  Mat<T> rs = Mat<T>::Zero(rows, cols);

  for (unsigned i = 0; i < tensor.size(); i++) {
    rs.row(i) = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>(
        tensor[i].data(), tensor[i].size());
  }
  return rs;
}

/*
 * Removed all the allocated arrays from the device
 */
template <typename T> EigenCuda<T>::~EigenCuda() {
  cublasDestroy(_handle);
  for (auto &p : _allocated)
    this->gpu_free(p.second);
}

/*
 * Allocate memory in the device using either pinned or pageable (default)
 * memory
 */
template <typename T> void EigenCuda<T>::gpu_alloc(T **x, std::size_t n) const {
  (_pinned) ? cudaMallocHost(x, n) : cudaMalloc(x, n);
}

/*
 * Deallocate memory from the device
 */
template <typename T> void EigenCuda<T>::gpu_free(T *x) const {
  (_pinned) ? cudaFreeHost(x) : cudaFree(x);
};

/*
 * Release the memory associated with the pointer `id` and removed the pointer
 * from the tracked pointers collection
 */
template <typename T> void EigenCuda<T>::free_matrix(int id) {
  // Free Array with id from the device
  gpu_free(_allocated.at(id));
  _allocated.erase(id);
}

/*
 * Allocate memory in the device for matrix A, then if if `copy_to_device`
 * copy the array to the device. Sometimes it only neccesary to allocate
 * space in the device without copying the array because the initial
 * values may not be important like a temporal matrix.
 */
template <typename T>
int EigenCuda<T>::initialize_Matrix(const Mat<T> &A, bool copy_to_device) {

  // size of the Matrices
  std::size_t size_A = A.size() * sizeof(T);

  // Pointer in the device
  T *dA;

  // Allocate either pageable or pinned memory
  gpu_alloc(&dA, size_A);

  // Track the allocated variables
  int id = _counter;
  _allocated.emplace(std::make_pair(id, dA));
  _counter += 1;

  // Transfer data to the GPU
  if (copy_to_device) {
    // Pointers at the host
    const T *hA = A.data();
    cudaMemcpy(dA, hA, size_A, cudaMemcpyHostToDevice);
  }

  return id;
}
/*
 * Call the gemm function from cublas, resulting in the multiplication of the
 * two matrices with identifiers id_A and id_B. The result is stored in
 * a Matrix (pointer) with identifier id_C.
 */
template <typename T>
void EigenCuda<T>::gemm(Shapes sh, std::tuple<int, int, int> ids) {
  // Invoke the gemm subroutine from cublas
  int id_A, id_B, id_C;
  std::tie(id_A, id_B, id_C) = ids;

  // Pointer to the arrays in the device
  T *dA, *dB, *dC;
  dA = _allocated.at(id_A);
  dB = _allocated.at(id_B);
  dC = _allocated.at(id_C);

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
 * Perform the matrix-matrix multiplication between A and B. First,
 * memory is allocated in the device for both matrices then a third temporal
 * array is allocated in the device that will contain the results. Finally, the
 * memory contains in the temporal result is copy back to the main memory and
 * Free the resources
 */
template <typename T>
Mat<T> EigenCuda<T>::dot(const Mat<T> &A, const Mat<T> &B) {
  // Matrix to store the result
  Mat<T> C = Mat<T>::Zero(A.rows(), B.cols());
  std::size_t size_C = C.size() * sizeof(T);

  // Indices of the matrices on the device
  std::tuple<int, int, int> ids = std::make_tuple(
      initialize_Matrix(A), initialize_Matrix(B), initialize_Matrix(C, false));

  // process on GPU
  Shapes sh{A.rows(), A.cols(), B.rows(), B.cols(), C.cols()};
  gemm(sh, ids);

  // send data back to CPU
  T *hC = C.data();
  T *dC = _allocated[std::get<2>(ids)];
  cudaMemcpy(hC, dC, size_C, cudaMemcpyDeviceToHost);

  // Free the result from the device
  free_matrix(std::get<0>(ids)); // Free A
  free_matrix(std::get<1>(ids)); // Free B
  free_matrix(std::get<2>(ids)); // Free C

  // create an eigen matrix
  C = Eigen::Map<Mat<T>>(hC, A.rows(), B.cols());

  return C;
}

template <typename T>
std::vector<Mat<T>>
EigenCuda<T>::triple_tensor_product(const Mat<T> &A, const Mat<T> &C,
                                    const std::vector<Mat<T>> &tensor) {
  // Perform the triple matrix multiplication A * matrix * C, for the vector
  // of matrices given by tensor
  std::vector<Mat<T>> rs(tensor.size());

  // Copy Matrix A and B to the device
  int id_A = initialize_Matrix(A);
  int id_C = initialize_Matrix(C);

  // allocate space in device for the temporal matrices
  int size_Y = A.rows() * C.cols() * sizeof(T);
  Mat<T> X = Mat<T>::Zero(A.cols(), C.cols());
  Mat<T> Y = Mat<T>::Zero(A.rows(), C.cols());
  Mat<T> matrix = Mat<T>::Zero(A.cols(), C.rows());

  int id_X = initialize_Matrix(X, false);
  int id_Y = initialize_Matrix(Y, false);
  int id_matrix = initialize_Matrix(matrix, false);

  // Iterate over the tensor Using the previous allocated space in the device
  transform(tensor.begin(), tensor.end(), rs.begin(),
            [this, id_A, id_C, id_X, id_Y, id_matrix, size_Y, &A, &C, &X,
             &Y](const Mat<T> &mtx) {
              assert(A.cols() == mtx.rows());
              assert(mtx.cols() == C.rows());

              // Copy matrix to the device
              T *d_matrix = _allocated.at(id_matrix);
              const T *h_mtx = mtx.data();

              // move temporal matrix to the preallocated space
              std::size_t size_mtx = mtx.size() * sizeof(T);
              cudaMemcpy(d_matrix, h_mtx, size_mtx, cudaMemcpyHostToDevice);

              // Compute first matrix multiplication
              Shapes sh1{mtx.rows(), mtx.cols(), C.rows(), C.cols(), X.rows()};
              std::tuple<int, int, int> ids =
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

  // Free all the allocated arrays from the device
  for (int x : {id_A, id_C, id_X, id_Y, id_matrix}) {
    free_matrix(x);
  }

  return rs;
}

/*
 * Multiply a matrix A by a 3D tensor represented as a vector of matrices.
 * Each iteration perform the operation mtx * A,  where mtx is the ith component
 * of the tensor
 */
template <typename T>
std::vector<Mat<T>>
EigenCuda<T>::right_matrix_tensor(const Mat<T> &A,
                                  const std::vector<Mat<T>> &tensor) {
  // result vector
  std::vector<Mat<T>> rs(tensor.size());

  // Copy Matrix A to the device
  int id_A = initialize_Matrix(A);

  // allocate space in device for the temporal matrices
  int rows = tensor[0].rows(); // rows of the submatrices
  int size_Y = rows * A.cols() * sizeof(T);
  Mat<T> Y = Mat<T>::Zero(rows, A.cols());
  Mat<T> matrix = Mat<T>::Zero(rows, A.rows());

  int id_Y = initialize_Matrix(Y, false);
  int id_matrix = initialize_Matrix(matrix, false);

  // Iterate over the tensor Using the previous allocated space in the device
  transform(tensor.begin(), tensor.end(), rs.begin(),
            [this, id_A, id_Y, id_matrix, size_Y, &A, &Y](const Mat<T> &mtx) {
              // Check that the matrix has the right dimension
              assert(mtx.cols() == A.rows());

              // Copy matrix to the device
              T *d_matrix = _allocated.at(id_matrix);
              const T *h_mtx = mtx.data();

              // move temporal matrix to the preallocated space
              std::size_t size_mtx = mtx.size() * sizeof(T);
              cudaMemcpy(d_matrix, h_mtx, size_mtx, cudaMemcpyHostToDevice);

              // Compute the matrix multiplication
              Shapes sh1{mtx.rows(), mtx.cols(), A.rows(), A.cols(), Y.rows()};
              std::tuple<int, int, int> ids =
                  std::make_tuple(id_matrix, id_A, id_Y);
              gemm(sh1, ids);

              // send data back to CPU
              T *hY = Y.data();
              T *dY = this->_allocated[id_Y];
              cudaMemcpy(hY, dY, size_Y, cudaMemcpyDeviceToHost);
              Y = Eigen::Map<Mat<T>>(hY, mtx.rows(), A.cols());

              return Y;
            });

  // Free all the allocated arrays from the device
  for (int x : {id_A, id_Y, id_matrix}) {
    free_matrix(x);
  }

  return rs;
}

// explicit instantiations
template class EigenCuda<float>;
template class EigenCuda<double>;
template Mat<float> stack<float>(const std::vector<Mat<float>> &);
template Mat<double> stack<double>(const std::vector<Mat<double>> &);
}  // namespace xtp
}  // namespace votca
