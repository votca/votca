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

  cpus_.resize(getNumberThreads());

#ifdef USE_CUDA
  Index no_gpus = count_available_gpus();
  gpus_.clear();
  if (inside_Parallel_region_) {
    if (threadID_parent_ < no_gpus) {
      gpus_.push_back(GPU_data(threadID_parent_));
    }
  } else {
    if (no_gpus > getNumberThreads()) {
      no_gpus = getNumberThreads();
    }
    for (Index i = 0; i < no_gpus; i++) {
      gpus_.push_back(GPU_data(i));
    }
  }
#endif
}

#ifdef USE_CUDA
void OpenMP_CUDA::setOperators(const std::vector<Eigen::MatrixXd>& tensor,
                               const Eigen::MatrixXd& rightoperator) {
  rOP_ = rightoperator;
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
  rOP_ = rightoperator;
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
    tensor *= rOP_();
  }
#else
  tensor *= rOP_();
#endif
  return;
}

void OpenMP_CUDA::setOperators(const Eigen::MatrixXd& leftoperator,
                               const Eigen::MatrixXd& rightoperator){
  lOP_ = leftoperator;
  rOP_ = rightoperator;
 
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
    matrix = lOP_() * matrix * rOP_();
  }
#else
  matrix = lOP_() * matrix * rOP_();
#endif
  return;
}
#ifdef USE_CUDA
void OpenMP_CUDA::createTemporaries(Index rows, Index cols) {

  std::for_each(cpus_.begin(), cpus_.end(),
                [&](CPU_data& d) { d.InitializeReduce(cols, cols); });

#pragma omp parallel for num_threads(gpus_.size())
  for (Index i = 0; i < Index(gpus_.size()); i++) {
    GPU_data& gpu = gpus_[i];
    gpu.activateGPU();
    gpu.push_back(rows, 1);
    gpu.push_back(rows, cols);
    gpu.push_back(rows, cols);
    gpu.push_back(cols, cols);
    gpu.temp.back()->setZero();
  }
}
#else
void OpenMP_CUDA::createTemporaries(Index, Index cols) {
  std::for_each(cpus_.begin(), cpus_.end(),
                [&](CPU_data& d) { d.InitializeReduce(cols, cols); });
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
    cpus_[threadid].ref_mat = matrix;
  }
#else
  cpus_[threadid].ref_mat = matrix;
#endif
}

void OpenMP_CUDA::A_TDA(const Eigen::VectorXd& vec) {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
  auto cpucomp = [&]() {
    CPU_data& cpu = cpus_[threadid];
    cpu.reduce() +=
        cpu.ref_mat().transpose() * vec.asDiagonal() * cpu.ref_mat();
  };
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.Mat(0).copy_to_gpu(vec);
    gpu.pipe().diag_gemm(gpu.Mat(1).transpose(), gpu.Mat(0), gpu.Mat(2));
    gpu.pipe().gemm(gpu.Mat(1).transpose(), gpu.Mat(2), gpu.Mat(3), 1.0);
  } else {
    cpucomp();
  }
#else
  cpucomp();
#endif
}

#ifdef USE_CUDA
void OpenMP_CUDA::createTemporaries(const Eigen::VectorXd& vec,
                                    const Eigen::MatrixXd& input, Index rows1,
                                    Index rows2, Index cols) {

  std::for_each(cpus_.begin(), cpus_.end(), [&](CPU_data& d) {
    d.InitializeReduce(input.rows(), input.cols());
    d.InitializeVec(input.rows());
  });

  rOP_ = input;
  vec_ = vec;

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
    gpu.push_back(reduction_[i].rows(), reduction_[i].cols());
    gpu.temp.back()->setZero();
  }
}

#else
void OpenMP_CUDA::createTemporaries(const Eigen::VectorXd& vec,
                                    const Eigen::MatrixXd& input, Index, Index,
                                    Index) {
  std::for_each(cpus_.begin(), cpus_.end(), [&](CPU_data& d) {
    d.InitializeReduce(input.rows(), input.cols());
    d.InitializeVec(input.rows());
  });

  rOP_ = input;
  vec_ = vec;
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
    cpus_[threadid].ref_mat = mat;
    mat *= vec_().asDiagonal();
  }
#else
  cpus_[threadid].ref_mat = mat;
  mat *= vec_().asDiagonal();
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
    cpus_[threadid].temp_vec.setZero();
  }
#else
  cpus_[threadid].temp_vec.setZero();
#endif
}

void OpenMP_CUDA::PrepareMatrix2(const Eigen::Block<const Eigen::MatrixXd>& mat,
                                 bool Hd2) {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
  auto cpucomp = [&]() {
    CPU_data& cpu = cpus_[threadid];
    if (Hd2) {
      Eigen::Map<Eigen::MatrixXd> row(cpu.temp_vec.data(), mat.rows(),
                                      cpu.ref_mat().rows());
      row += mat * cpus_[threadid].ref_mat().transpose();
    } else {
      Eigen::Map<Eigen::MatrixXd> row(cpu.temp_vec.data(), cpu.ref_mat().rows(),
                                      mat.rows());
      row += cpu.ref_mat() * mat.transpose();
    }
  };
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
    cpucomp();
  }
#else
  cpucomp();
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
    cpus_[threadid].temp_vec += row;
  }
#else
  cpus_[threadid].temp_vec += row;
#endif
}

void OpenMP_CUDA::MultiplyRow(Index row) {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
  auto cpucomp = [&]() {
    cpus_[threadid].reduce().row(row) =
        cpus_[threadid].temp_vec.transpose() * rOP_();
  };
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.pipe().gemm(gpu.Mat(4).transpose(), gpu.Mat(1),
                    gpu.Mat(6).block(row, 0, 1, gpu.Mat(6).cols()), 0.0);
  } else {
    cpucomp();
  }
#else
  cpucomp();
#endif
}

#ifdef USE_CUDA
void OpenMP_CUDA::createAdditionalTemporaries(Index rows, Index cols) {
#pragma omp parallel for num_threads(gpus_.size())
  for (Index i = 0; i < Index(gpus_.size()); i++) {
    GPU_data& gpu = gpus_[i];
    gpu.activateGPU();
    gpu.resize(2, rows, cols);
    gpu.resize(3, rows, cols);
    gpu.resize(4, rows, rows);
  }
}
#else
void OpenMP_CUDA::createAdditionalTemporaries(Index, Index) { ; }
#endif

void OpenMP_CUDA::PushMatrix1(Eigen::MatrixXd& mat) {

  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.Mat(2).copy_to_gpu(mat);
  } else {
    cpus_[threadid].ref_mat = mat;
  }
#else
  cpus_[threadid].ref_mat = mat;
#endif
}

void OpenMP_CUDA::MultiplyBlocks(const Eigen::Block<const Eigen::MatrixXd>& mat,
                                 Index i1, Index i2) {
  Index parentid = getParentThreadId();
  Index threadid = getLocalThreadId(parentid);
  auto cpucomp = [&]() {
    CPU_data& cpu = cpus_[threadid];
    Eigen::MatrixXd block = cpu.ref_mat() * mat.transpose();
    cpu.reduce().block(i1 * block.rows(), 0, block.rows(),
                       cpu.reduce().cols()) +=
        block * rOP_().block(i2 * block.rows(), 0, block.rows(), rOP_().cols());
    if (i1 != i2) {
      cpu.reduce().block(i2 * block.rows(), 0, block.rows(),
                         cpu.reduce().cols()) +=
          block.transpose() *
          rOP_().block(i1 * block.rows(), 0, block.rows(), rOP_().cols());
    }
  };
#ifdef USE_CUDA
  if (isGPUthread(parentid)) {
    GPU_data& gpu = gpus_[threadid];
    gpu.activateGPU();
    gpu.Mat(3).copy_to_gpu(mat);
    gpu.pipe().gemm(gpu.Mat(2), gpu.Mat(3).transpose(), gpu.Mat(4));
    Index blocksize = gpu.Mat(4).rows();
    Index inputcols = gpu.Mat(1).cols();
    gpu.pipe().gemm(
        gpu.Mat(4), gpu.Mat(1).block(i2 * blocksize, 0, blocksize, inputcols),
        gpu.Mat(6).block(i1 * blocksize, 0, blocksize, inputcols), 1.0);
    if (i1 != i2) {
      gpu.pipe().gemm(gpu.Mat(4).transpose(),
                      gpu.Mat(1).block(i1 * blocksize, 0, blocksize, inputcols),
                      gpu.Mat(6).block(i2 * blocksize, 0, blocksize, inputcols),
                      1.0);
    }
  } else {
    cpucomp();
  }
#else
  cpucomp();
#endif
}

Eigen::MatrixXd OpenMP_CUDA::getReductionVar() {
#ifdef USE_CUDA
#pragma omp parallel for num_threads(gpus_.size())
  for (Index i = 0; i < Index(gpus_.size()); i++) {
    GPU_data& gpu = gpus_[i];
    gpu.activateGPU();
    cpus_[i].reduce() = *(gpu.temp.back());
  }
#endif
  for (Index i = 1; i < Index(cpus_.size()); i++) {
    cpus_[0].reduce() += cpus_[i].reduce();
  }
  return cpus_[0].reduce();
}

}  // namespace xtp
}  // namespace votca
