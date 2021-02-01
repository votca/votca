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

  inside_Parallel_region_ = OPENMP::InsideActiveParallelRegion();
  threadID_parent_ = OPENMP::getThreadId();

  cpu_intermediate_input_.resize(getNumberThreads());
  cpu_temporaries_.resize(getNumberThreads());
#ifdef USE_CUDA
  Index no_gpus = count_available_gpus();
  gpus_.clear();
  if (inside_Parallel_region_) {
    if (threadID_parent_ < no_gpus) {
      gpus_.push_back(GPU_data(threadID_parent_));
    }
  } else {
    for (Index i = 0; i < no_gpus; i++) {
      gpus_.push_back(GPU_data(i));
    }
  }
#endif
}

#ifdef USE_CUDA
void OpenMP_CUDA::setOperators(const std::vector<Eigen::MatrixXd>& tensor,
                               const Eigen::MatrixXd& rightoperator) {
  rightoperator_ = &rightoperator;
#pragma omp parallel for num_threads(gpus_.size())
  for (Index i = 0; i < Index(gpus_.size()); i++) {
    GPU_data& gpu = gpus_[i];
    gpu.activateGPU();
    const Eigen::MatrixXd& head = tensor.front();
    gpu.push_back(head.rows(), head.cols());
    gpu.push_back(rightoperator);
    gpu.push_back(head.rows(), rightoperator.cols());
  }
}
#else
void OpenMP_CUDA::setOperators(const std::vector<Eigen::MatrixXd>&,
                               const Eigen::MatrixXd& rightoperator) {
  rightoperator_ = &rightoperator;
}
#endif

#ifdef USE_CUDA

bool OpenMP_CUDA::isInVector(Index Id, const std::vector<GPU_data>& vec) {
  return (std::find_if(vec.begin(), vec.end(), [&Id](const GPU_data& d) {
            return d.Id == Id;
          }) != vec.end());
}
bool OpenMP_CUDA::isGPUthread(Index ParentThreadId) const {
  return isInVector(ParentThreadId, gpus_);
}
#endif

Index OpenMP_CUDA::getParentThreadId() const {
  return inside_Parallel_region_ ? threadID_parent_ : OPENMP::getThreadId();
}
Index OpenMP_CUDA::getLocalThreadId(Index ParentThreadId) const {
  return inside_Parallel_region_ ? 0 : ParentThreadId;
}
Index OpenMP_CUDA::getNumberThreads() const {
  return inside_Parallel_region_ ? 1 : OPENMP::getMaxThreads();
}

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
  Index threadid = getParentThreadId();
  if (isGPUthread(threadid)) {
    GPU_data& gpu = gpus_[getLocalThreadId(threadid)];
    gpu.activateGPU();
    gpu.Mat(0).copy_to_gpu(tensor);
    gpu.pipe().gemm(gpu.Mat(0), gpu.Mat(1), gpu.Mat(2));
    tensor = gpu.Mat(2);
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
#pragma omp parallel for num_threads(gpus_.size())
  for (Index i = 0; i < Index(gpus_.size()); i++) {
    GPU_data& gpu = gpus_[i];
    gpu.activateGPU();
    gpu.push_back(leftoperator);
    gpu.push_back(leftoperator.cols(), rightoperator.rows());
    gpu.push_back(leftoperator.rows(), rightoperator.rows());
    gpu.push_back(rightoperator);
    gpu.push_back(leftoperator.rows(), rightoperator.cols());
  }
#endif
}

void OpenMP_CUDA::MultiplyLeftRight(Eigen::MatrixXd& matrix) {
#ifdef USE_CUDA
  Index threadid = getParentThreadId();
  if (isGPUthread(threadid)) {
    GPU_data& gpu = gpus_[getLocalThreadId(threadid)];
    gpu.activateGPU();
    gpu.Mat(1).copy_to_gpu(matrix);
    gpu.pipe().gemm(gpu.Mat(0), gpu.Mat(1), gpu.Mat(2));
    gpu.pipe().gemm(gpu.Mat(2), gpu.Mat(3), gpu.Mat(4));
    matrix = gpu.Mat(4);
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
  reduction_ = std::vector<Eigen::MatrixXd>(getNumberThreads(),
                                            Eigen::MatrixXd::Zero(cols, cols));

#pragma omp parallel for num_threads(gpus_.size())
  for (Index i = 0; i < Index(gpus_.size()); i++) {
    GPU_data& gpu = gpus_[i];
    gpu.activateGPU();
    gpu.push_back(rows, 1);
    gpu.push_back(rows, cols);
    gpu.push_back(rows, cols);
    gpu.push_back(reduction_[i]);
  }
}
#else
void OpenMP_CUDA::createTemporaries(Index, Index cols) {
  reduction_ = std::vector<Eigen::MatrixXd>(getNumberThreads(),
                                            Eigen::MatrixXd::Zero(cols, cols));
}

#endif

void OpenMP_CUDA::PushMatrix(Eigen::MatrixXd& matrix) {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.Mat(1).copy_to_gpu(matrix);
  } else {
    cpu_intermediate_input_[threadid] = &matrix;
  }
#else
  cpu_intermediate_input_[threadid] = &matrix;
#endif
}

void OpenMP_CUDA::A_TDA(const Eigen::VectorXd& vec) {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.Mat(0).copy_to_gpu(vec);
    gpu.pipe().diag_gemm(gpu.Mat(1), gpu.Mat(0), gpu.Mat(2));
    gpu.pipe().gemm(gpu.Mat(1).transpose(), gpu.Mat(2), gpu.Mat(3), 1.0);
  } else {
    reduction_[threadid] += cpu_intermediate_input_[threadid]->transpose() *
                            vec.asDiagonal() *
                            (*cpu_intermediate_input_[threadid]);
  }
#else
  reduction_[threadid] += cpu_intermediate_input_[threadid]->transpose() *
                          vec.asDiagonal() *
                          (*cpu_intermediate_input_[threadid]);
#endif
}

#ifdef USE_CUDA
void OpenMP_CUDA::createTemporaries(const Eigen::VectorXd& vec,
                                    const Eigen::MatrixXd& input, Index rows1,
                                    Index rows2, Index cols) {
  reduction_ = std::vector<Eigen::MatrixXd>(
      getNumberThreads(), Eigen::MatrixXd::Zero(input.rows(), input.cols()));
  temp_ = std::vector<Eigen::VectorXd>(
      getNumberThreads(), Eigen::VectorXd::Zero(input.cols()));
  rightoperator_ = &input;
  vec_ = &vec;

#pragma omp parallel for num_threads(gpus_.size())
  for (Index i = 0; i < Index(gpus_.size()); i++) {
    GPU_data& gpu = gpus_[i];
    gpu.activateGPU();
    gpu.push_back(vec);
    gpu.push_back(input);
    gpu.push_back(rows1, cols);
    gpu.push_back(rows2, cols);
    gpu.push_back(rows1 * rows2, 1);
    gpu.push_back(rows1 * rows2, 1);
    gpu.push_back(reduction_[i]);
  }
}

#else
void OpenMP_CUDA::createTemporaries(const Eigen::VectorXd& vec,
                                    const Eigen::MatrixXd& input, Index, Index,
                                    Index) {
  reduction_ = std::vector<Eigen::MatrixXd>(
      getNumberThreads(), Eigen::MatrixXd::Zero(input.rows(), input.cols()));
  temp_ = std::vector<Eigen::VectorXd>(
      getNumberThreads(), Eigen::VectorXd::Zero(input.rows()));
  rightoperator_ = &input;
  vec_ = &vec;
}
#endif

void OpenMP_CUDA::PrepareMatrix1(Eigen::MatrixXd& mat) {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.Mat(2).copy_to_gpu(mat);
    gpu.pipe().diag_gemm(gpu.Mat(2), gpu.Mat(0), gpu.Mat(2));
  } else {
    cpu_intermediate_input_[threadid] = &mat;
    mat *= vec_->asDiagonal();
  }
#else
  cpu_intermediate_input_[threadid] = &mat;
  mat *= vec_->asDiagonal();
#endif
}

void OpenMP_CUDA::SetTempZero() {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.Mat(4).setZero();
  } else {
    temp_[threadid].setZero();
  }
#else
  temp_[threadid].setZero();
#endif
}

void OpenMP_CUDA::PrepareMatrix2(const Eigen::MatrixXd& mat, bool Hd2) {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.Mat(3).copy_to_gpu(mat);
    if (Hd2) {
      gpu.Mat(4).reshape(gpu.Mat(3).rows(), gpu.Mat(2).rows());
      gpu.pipe().gemm(gpu.Mat(3), gpu.Mat(2).transpose(), gpu.Mat(4), 1.0);
    } else {
      gpu.Mat(4).reshape(gpu.Mat(2).rows(), gpu.Mat(3).rows());
      gpu.pipe().gemm(gpu.Mat(2), gpu.Mat(3).transpose(), gpu.Mat(4), 1.0);
    }
    gpu.Mat(4).reshape(gpu.Mat(2).rows() * gpu.Mat(3).rows(), 1);
  } else {

    if (Hd2) {
      Eigen::Map<Eigen::MatrixXd> row(
          temp_[threadid].data(), mat.rows(),
          cpu_intermediate_input_[threadid]->rows());
      row += mat * cpu_intermediate_input_[threadid]->transpose();
    } else {
      Eigen::Map<Eigen::MatrixXd> row(temp_[threadid].data(),
                                      cpu_intermediate_input_[threadid]->rows(),
                                      mat.rows());
      row += *cpu_intermediate_input_[threadid] * mat.transpose();
    }
  }
#else
  if (Hd2) {
    Eigen::Map<Eigen::MatrixXd> row(temp_[threadid].data(), mat.rows(),
                                    cpu_intermediate_input_[threadid]->rows());
    row += mat * cpu_intermediate_input_[threadid]->transpose();
  } else {
    Eigen::Map<Eigen::MatrixXd> row(temp_[threadid].data(),
                                    cpu_intermediate_input_[threadid]->rows(),
                                    mat.rows());
    row += *cpu_intermediate_input_[threadid] * mat.transpose();
  }
#endif
}

void OpenMP_CUDA::Addvec(const Eigen::VectorXd& row) {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.Mat(5).copy_to_gpu(row);
    gpu.pipe().axpy(gpu.Mat(5), gpu.Mat(4), 1.0);
  } else {
    temp_[threadid] += row;
  }
#else
  temp_[threadid] += row;
#endif
}

void OpenMP_CUDA::MultiplyRow(Index row) {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.pipe().gemm(gpu.Mat(4).transpose(), gpu.Mat(1),
                    gpu.Mat(6).block(row, 0, 1, gpu.Mat(1).cols()), 0.0);
  } else {
    reduction_[threadid].row(row) = temp_[threadid].transpose() * (*rightoperator_);
  }
#else
  reduction_[threadid].row(row) = temp_[threadid].transpose() * (*rightoperator_);
#endif
}

Eigen::MatrixXd OpenMP_CUDA::getReductionVar() {
#ifdef USE_CUDA
#pragma omp parallel for num_threads(gpus_.size())
  for (Index i = 0; i < Index(gpus_.size()); i++) {
    GPU_data& gpu = gpus_[i];
    gpu.activateGPU();
    reduction_[i] = *(gpu.temp.back());
  }
#endif
  for (Index i = 1; i < Index(reduction_.size()); i++) {
    reduction_[0] += reduction_[i];
  }
  return reduction_[0];
}

}  // namespace xtp
}  // namespace votca
