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
 * Stack a vector of matrices as a single matrix, where each column corresponds
 * to a matrix.
 */
template <typename T>
Mat<T> stack(const std::vector<Mat<T>> &tensor) {

  int rows = tensor[0].size();  // size of each matrix
  int cols = tensor.size();     // number of matrices in tensor

  // row major to save the tensor
  Mat<T> rs = Mat<T>::Zero(rows, cols);

  for (auto i = 0; i < cols; i++) {
    rs.col(i) = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>(
        tensor[i].data(), tensor[i].size());
  }
  return rs;
}

/*
 * Removed all the allocated arrays from the device
 */
template <typename T>
EigenCuda<T>::~EigenCuda() {
  cublasDestroy(_handle);
  for (auto &p : _allocated) this->gpu_free(p.second);
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
  // Deallocate memory from the device
  (_pinned) ? cudaFreeHost(x) : cudaFree(x);
};

/*
 * Release the memory associated with the pointer `id` and removed the pointer
 * from the tracked pointers collection
 */
template <typename T>
void EigenCuda<T>::free_matrix(int id) {
  // Free Array with id from the device
  gpu_free(_allocated.at(id));
  _allocated.erase(id);
}

/*
 * Method to shift pointers of the allocated tensor in the device. When
 * iterating the tensor this method is invoked to get the next submatrix from a
 * given tensor
 */
template <typename T>
void EigenCuda<T>::shift_pointers_by(const std::vector<int> &pointers,
                                     const std::vector<long int> &shifts) {
  for (unsigned i = 0; i < pointers.size(); i++) {
    _allocated.at(pointers[i]) += static_cast<int>(shifts[i]);
  }
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
 * see: https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-gemm
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
 * the resources are free
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
  free_matrix(std::get<0>(ids));  // Free A
  free_matrix(std::get<1>(ids));  // Free B
  free_matrix(std::get<2>(ids));  // Free C

  // create an eigen matrix
  C = Eigen::Map<Mat<T>>(hC, A.rows(), B.cols());

  return C;
}

/*
 * \brief performs a matrix_1 * tensor * matrix_2 multiplication
 * \return matrix where each column is the result of the matrices
 * multiplication.
 *
 * Initially, it allocates memory and copy the matrices A and C together with
 * the tensor to the device. Also, the function allocates the result tensor Y
 * and a temporal matrix X.
 * This last matrix is not copy into the device because is initial value is not
 * relevant. Subsequently, the method iterates over each submatrix in `tensor`
 * and perform the following operations: X = tensor(i) * C Y(i) = A * X then the
 * final Y is copy back to main memory. This final matrix Y contains in each
 * column the result of the tensor operation. Also, notice that the matrix X is
 * never set to zero after each iteration because the gemm function perform the
 * matrix multiplication: R = alpha M * N + beta R where alpha and beta are two
 * scalar constants set to 1 and 0 respectively. Therefore, X is ALWAYS SET TO
 * ZERO BEFORE THE MATRIX MULTIPLICATION.
 */

template <typename T>
Mat<T> EigenCuda<T>::triple_tensor_product(const Mat<T> &A, const Mat<T> &C,
                                           const std::vector<Mat<T>> &tensor) {
  // Copy Matrix A and B to the device
  int id_A = initialize_Matrix(A);
  int id_C = initialize_Matrix(C);

  // Stack tensor into a single matrix
  Mat<T> super_matrix = eigencuda::stack(tensor);

  // Move tensor to the device
  int id_super = initialize_Matrix(super_matrix);

  // rows and cols of matrices store in the tensor
  int mtx_rows = tensor[0].rows();
  int mtx_cols = tensor[0].cols();

  // Allocate space fo the result Tensor
  Mat<T> Y = Mat<T>::Zero(mtx_rows * C.cols(), super_matrix.cols());
  int id_Y = initialize_Matrix(Y, false);
  T *init_Y = _allocated.at(id_Y);  // pointer to initial location

  // allocate space in device for the temporal matrices
  Mat<T> X = Mat<T>::Zero(A.cols(), C.cols());
  int id_X = initialize_Matrix(X, false);

  // Iterate over the tensor Using the previous allocated space in the device
  for (unsigned i = 0; i < tensor.size(); i++) {

    // Compute first matrix multiplication
    Shapes sh1{mtx_rows, mtx_cols, C.rows(), C.cols(), X.rows()};
    std::tuple<int, int, int> ids = std::make_tuple(id_super, id_C, id_X);
    gemm(sh1, ids);

    // compute the second matrix multiplication
    Shapes sh2{A.rows(), A.cols(), X.rows(), X.cols(), C.cols()};
    ids = std::make_tuple(id_A, id_X, id_Y);
    gemm(sh2, ids);

    // shift the pointer containing the super_matrix and the result tensor
    std::vector<int> pointers{id_super, id_Y};
    std::vector<long int> shifts{mtx_rows * mtx_cols, mtx_rows * C.cols()};
    shift_pointers_by(pointers, shifts);
  }

  // send data back to CPU
  T *hY = Y.data();
  size_t size_Y = Y.size() * sizeof(T);
  cudaMemcpy(hY, init_Y, size_Y, cudaMemcpyDeviceToHost);
  Y = Eigen::Map<Mat<T>>(hY, Y.rows(), Y.cols());

  // Free all the allocated arrays from the device
  for (int x : {id_A, id_C, id_X, id_Y, id_super}) {
    free_matrix(x);
  }

  return Y;  // each column contains the resulting product
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
Mat<T> EigenCuda<T>::right_matrix_tensor(const Mat<T> &A,
                                         const std::vector<Mat<T>> &tensor) {
  // Copy Matrix A to the device
  int id_A = initialize_Matrix(A);

  // Stack tensor into a single matrix
  Mat<T> super_matrix = eigencuda::stack(tensor);

  // Move tensor to the device
  int id_super = initialize_Matrix(super_matrix);

  // rows and cols of matrices store in the tensor
  int mtx_rows = tensor[0].rows();
  int mtx_cols = tensor[0].cols();

  // Allocate space fo the result Tensor
  Mat<T> Y = Mat<T>::Zero(mtx_rows * A.cols(), super_matrix.cols());
  int id_Y = initialize_Matrix(Y, false);
  T *init_Y = _allocated.at(id_Y);  // pointer to initial location

  // Iterate over the tensor Using the previous allocated space in the device
  for (unsigned i = 0; i < tensor.size(); i++) {

    // Compute the matrix multiplication
    Shapes sh1{mtx_rows, mtx_cols, A.rows(), A.cols(), mtx_rows};
    std::tuple<int, int, int> ids = std::make_tuple(id_super, id_A, id_Y);
    gemm(sh1, ids);

    // shift the pointer containing the super_matrix and the result tensor
    std::vector<int> pointers{id_super, id_Y};
    std::vector<long int> shifts{mtx_rows * mtx_cols, mtx_rows * A.cols()};
    shift_pointers_by(pointers, shifts);
  }

  // send data back to CPU
  T *hY = Y.data();
  size_t size_Y = Y.size() * sizeof(T);
  cudaMemcpy(hY, init_Y, size_Y, cudaMemcpyDeviceToHost);
  Y = Eigen::Map<Mat<T>>(hY, Y.rows(), Y.cols());

  // Free all the allocated arrays from the device
  for (int x : {id_A, id_Y, id_super}) {
    free_matrix(x);
  }

  return Y;
}

// explicit instantiations
template class EigenCuda<float>;
template class EigenCuda<double>;
template Mat<float> stack<float>(const std::vector<Mat<float>> &);
template Mat<double> stack<double>(const std::vector<Mat<double>> &);

}  // namespace xtp
}  // namespace votca
