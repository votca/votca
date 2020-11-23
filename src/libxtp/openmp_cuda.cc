/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Local VOTCA includes
#include "votca/xtp/openmp_cuda.h"

namespace votca {
namespace xtp {

OpenMP_CUDA::OpenMP_CUDA() {
#ifdef USE_CUDA
  no_gpus_ = count_available_gpus();
  for (Index i = 0; i < no_gpus_; i++) {
    cuda_pips_.push_back(std::make_unique<CudaPipeline>(int(i)));
  }
  temp_ = std::vector<temporaries>(no_gpus_);
#endif
}

#ifdef USE_CUDA
void OpenMP_CUDA::setOperators(const std::vector<Eigen::MatrixXd>& tensor,
                               const Eigen::MatrixXd& rightoperator) {
  rightoperator_ = &rightoperator;
#pragma omp parallel for num_threads(no_gpus_)
  for (Index i = 0; i < no_gpus_; i++) {
    const Eigen::MatrixXd& head = tensor.front();
    const cudaStream_t& stream = cuda_pips_[i]->get_stream();
    temp_[i].A = std::make_unique<CudaMatrix>(head.rows(), head.cols(), stream);
    temp_[i].B = std::make_unique<CudaMatrix>(rightoperator, stream);
    temp_[i].C =
        std::make_unique<CudaMatrix>(head.rows(), rightoperator.cols(), stream);
  }
}
#else
void OpenMP_CUDA::setOperators(const std::vector<Eigen::MatrixXd>&,
                               const Eigen::MatrixXd& rightoperator) {
  rightoperator_ = &rightoperator;
}
#endif

/*
 * The Cuda device behaves like a server that is receiving matrix-matrix
 * multiplications from a single stream (an Nvidia queue) and handle them
 * in an asynchronous way. It performs the following operations when recieving
 * a request:
 *  1. Check that there is enough space for the arrays
 *  2. Allocate memory for each matrix
 *  3. Copy the matrix to the allocated space
 *  4. Perform the matrix multiplication
 *  5. Return the result matrix
 * The Cuda device knows to which memory address it needs to copy back the
 * result. see: https://docs.nvidia.com/cuda/cublas/index.html#thread-safety2
 */
void OpenMP_CUDA::MultiplyRight(Eigen::MatrixXd& tensor) {

  // All the GPU communication happens through a single thread that reuses
  // all memory allocated in the GPU and it's dynamically load-balanced by
  // OpenMP. The rest of the threads use the default CPU matrix
  // multiplication
#ifdef USE_CUDA
  Index threadid = OPENMP::getThreadId();
  if (threadid < no_gpus_) {
    temp_[threadid].A->copy_to_gpu(tensor);
    cuda_pips_[threadid]->gemm(*temp_[threadid].A, *temp_[threadid].B,
                              *temp_[threadid].C);
    tensor = *temp_[threadid].C;
  } else {
    tensor *= *rightoperator_;
  }
#else
  tensor *= *rightoperator_;
#endif
  return;
}

void OpenMP_CUDA::setOperators(const Eigen::MatrixXd& leftoperator,
                               const Eigen::MatrixXd& rightoperator) {
  rightoperator_ = &rightoperator;
  leftoperator_ = &leftoperator;
#ifdef USE_CUDA
#pragma omp parallel for num_threads(no_gpus_)
  for (Index i = 0; i < no_gpus_; i++) {
    const cudaStream_t& stream = cuda_pips_[i]->get_stream();
    temp_[i].A = std::make_unique<CudaMatrix>(leftoperator, stream);
    temp_[i].B = std::make_unique<CudaMatrix>(leftoperator.cols(),
                                              rightoperator.rows(), stream);
    temp_[i].C = std::make_unique<CudaMatrix>(leftoperator.rows(),
                                              rightoperator.rows(), stream);
    temp_[i].D = std::make_unique<CudaMatrix>(rightoperator, stream);
    temp_[i].E = std::make_unique<CudaMatrix>(leftoperator.rows(),
                                              rightoperator.cols(), stream);
  }
#endif
}

void OpenMP_CUDA::MultiplyLeftRight(Eigen::MatrixXd& matrix) {
#ifdef USE_CUDA
  Index threadid = OPENMP::getThreadId();
  if (threadid < no_gpus_) {
    temp_[threadid].B->copy_to_gpu(matrix);
    cuda_pips_[threadid]->gemm(*temp_[threadid].A, *temp_[threadid].B,
                              *temp_[threadid].C);
    cuda_pips_[threadid]->gemm(*temp_[threadid].C, *temp_[threadid].D,
                              *temp_[threadid].E);
    matrix = *temp_[threadid].E;
  } else {
    matrix = (*leftoperator_) * matrix * (*rightoperator_);
  }
#else
  matrix = (*leftoperator_) * matrix * (*rightoperator_);
#endif
  return;
}
#ifdef USE_CUDA
void OpenMP_CUDA::createTemporaries(Index rows, Index cols) {
  reduction_ = std::vector<Eigen::MatrixXd>(OPENMP::getMaxThreads(),
                                            Eigen::MatrixXd::Zero(cols, cols));

#pragma omp parallel for num_threads(no_gpus_)
  for (Index i = 0; i < no_gpus_; i++) {
    const cudaStream_t& stream = cuda_pips_[i]->get_stream();
    temp_[i].A = std::make_unique<CudaMatrix>(rows, 1, stream);
    temp_[i].B = std::make_unique<CudaMatrix>(rows, cols, stream);
    temp_[i].C = std::make_unique<CudaMatrix>(rows, cols, stream);
    temp_[i].D = std::make_unique<CudaMatrix>(reduction_[0], stream);
  }
}
#else
void OpenMP_CUDA::createTemporaries(Index, Index cols) {
  reduction_ = std::vector<Eigen::MatrixXd>(OPENMP::getMaxThreads(),
                                            Eigen::MatrixXd::Zero(cols, cols));
}

#endif

void OpenMP_CUDA::A_TDA(const Eigen::MatrixXd& matrix,
                        const Eigen::VectorXd& vec) {
  Index threadid = OPENMP::getThreadId();
#ifdef USE_CUDA
  if (threadid < no_gpus_) {
    temp_[threadid].A->copy_to_gpu(vec);
    temp_[threadid].B->copy_to_gpu(matrix);
    cuda_pips_[threadid]->diag_gemm(*temp_[threadid].B, *temp_[threadid].A,
                                   *temp_[threadid].C);
    cuda_pips_[threadid]->gemm(*temp_[threadid].B, *temp_[threadid].C,
                              *temp_[threadid].D, true, false, 1.0);
  } else {
    reduction_[threadid] += matrix.transpose() * vec.asDiagonal() * matrix;
  }
#else
  reduction_[threadid] += matrix.transpose() * vec.asDiagonal() * matrix;
#endif
}

Eigen::MatrixXd OpenMP_CUDA::A_TDA_result() {
#ifdef USE_CUDA
#pragma omp parallel for num_threads(no_gpus_)
  for (Index i = 0; i < no_gpus_; i++) {
    reduction_[i] = *temp_[i].D;
  }
#endif
  for (Index i = 1; i < Index(reduction_.size()); i++) {
    reduction_[0] += reduction_[i];
  }
  return reduction_[0];
}

}  // namespace xtp
}  // namespace votca
