/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __XTP_EIGEN_CUDA__H
#define __XTP_EIGEN_CUDA__H

#include <cublas_v2.h>
#include <curand.h>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <votca/xtp/eigen.h>

/**
 * \brief Perform matrix-matrix multiplication in a GPU
 *
 * The `EigenCuda` class handles the allocation and deallocation of arrays on
 * the GPU. Firstly, to perform a matrix multiplication, memory must be
 * allocated in the device to contain the involved matrices. The
 * `initialize_matrix` method firstly allocates memory by calling the
 * `gpu_alloc` method that allocates either pinned or pageable memory, see:
 * https://devblogs.nvidia.com/how-optimize-data-transfers-cuda-cc/ Then the
 * array could be optionally copy to the device. The `initialize_matrix` method
 * internally store the pointer in the device to the allocated array and returns
 * and identifier that represents such pointer. This identifier/pointer
 * (key/value) mechanism allows to reuse already allocated space in the device
 * by bookkeeping the identifiers representing the memory.
 *
 */

namespace votca {
namespace xtp {

// Structure with the sizes to call ?GEMM
struct Shapes {
  int A_rows;
  int A_cols;
  int B_rows;
  int B_cols;
  int C_rows;

  Shapes(long int _a_rows, long int _a_cols, long int _b_rows, long int _b_cols,
         long int _c_rows)
      : A_rows{static_cast<int>(_a_rows)},
        A_cols{static_cast<int>(_a_cols)},
        B_rows{static_cast<int>(_b_rows)},
        B_cols{static_cast<int>(_b_cols)},
        C_rows{static_cast<int>(_c_rows)} {}
};
// col Major for CUDA
template <typename T>
using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

// Class to manage the GPU resourcess to perform matrix-matrix multiplications
template <typename T>
class EigenCuda {

 public:
  EigenCuda() { cublasCreate(&_handle); }
  EigenCuda(bool pinned) : _pinned{pinned} { cublasCreate(&_handle); }

  // Deallocate both the handler and allocated arrays
  ~EigenCuda();

  // Remove the copy operations
  EigenCuda(const EigenCuda &) = delete;
  EigenCuda &operator=(const EigenCuda &) = delete;

  // Matrix matrix multiplication
  Mat<T> dot(const Mat<T> &A, const Mat<T> &B);

  // Perform the triple matrix multiplication A * matrix * C, for the vector
  // of matrices given by tensor
  std::vector<Mat<T>> triple_tensor_product(const Mat<T> &A, const Mat<T> &C,
                                            const std::vector<Mat<T>> &tensor);

  // Perform a multiplication between a matrix and a tensor
  std::vector<Mat<T>> right_matrix_tensor(const Mat<T> &A,
                                          const std::vector<Mat<T>> &tensor);

 private:
  // Allocate memory in the device
  void gpu_alloc(T **x, std::size_t n) const;

  // Deallocate memory from the device
  void gpu_free(T *x) const;

  // Allocate memory in the device, optionally copying the array to the GPU
  int initialize_Matrix(const Mat<T> &A, bool copy_to_device = true);

  // Invoke the ?gemm function of cublas
  void gemm(Shapes shapes, std::tuple<int, int, int> ids);

  // Deallocate Matrix identifier `id` from the device
  void free_matrix(int id);

  // Cuda variables
  cublasHandle_t _handle;
  bool _pinned = false;

  // Allocation booking
  int _counter = 0;
  std::unordered_map<int, T *> _allocated;
};

// Stack a vector of matrices as a matrix where is row contains a matrix
template <typename T>
Mat<T> stack(const std::vector<Mat<T>> &tensor);

}  // namespace xtp
}  // namespace votca

#endif /*XTP_EIGEN_CUDA__H */
