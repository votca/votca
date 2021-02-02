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

#pragma once
#ifndef VOTCA_XTP_OPENMP_CUDA_H
#define VOTCA_XTP_OPENMP_CUDA_H

// Local VOTCA includes
#include "eigen.h"

#ifdef USE_CUDA
#include "cudapipeline.h"
#endif

/**
 * \brief Supports operations on Matrices using OPENMP and
 * CUDA.
 *
 * Each operation works with 2-3 steps
 * 1) Allocate temporary matrices and move fixed data to the gpu before the
 * openmp region is created 2) Inside the openmp region, move the loop data to
 * the GPU and perform calculation there 3) For reduction operations, transfer
 * the GPU data back to the CPU after the loop is finished Each GPU is served by
 * one CPU thread, the other CPU threads perform the normal CPU based operations
 * If no GPU is present all CPUs simply do CPU work.
 * If this class is created inside an OPENMP region, it still ensures, that over
 * that OPENMP region not more threads access the GPUs then GPUs are present.
 * Otherwise it will work purely in serial. So this class does NOT work with
 * nested OPENMP
 */

namespace votca {
namespace xtp {

class OpenMP_CUDA {
 public:
  OpenMP_CUDA();
  static Index UsingGPUs() {
#ifdef USE_CUDA
    return count_available_gpus();
#else
    return 0;
#endif
  }

  // 3c multiply
  void setOperators(const std::vector<Eigen::MatrixXd>& tensor,
                    const Eigen::MatrixXd& rightoperator);
  void MultiplyRight(Eigen::MatrixXd& matrix);

  // 3c
  void setOperators(const Eigen::MatrixXd& leftoperator,
                    const Eigen::MatrixXd& rightoperator);
  void MultiplyLeftRight(Eigen::MatrixXd& matrix);

  // RPA
  void createTemporaries(Index rows, Index cols);
  void PushMatrix(Eigen::MatrixXd& mat);
  void A_TDA(const Eigen::VectorXd& vec);

  // Hd + Hqp + Hd2
  void createTemporaries(const Eigen::VectorXd& vec,
                         const Eigen::MatrixXd& input, Index rows1, Index rows2,
                         Index cols);
  void PrepareMatrix1(Eigen::MatrixXd& mat);
  void SetTempZero();
  void PrepareMatrix2(const Eigen::Block<const Eigen::MatrixXd>& mat, bool Hd2);
  void Addvec(const Eigen::VectorXd& row);
  void MultiplyRow(Index row);

  // Hx

  void createAdditionalTemporaries(Index rows, Index cols);
  void PushMatrix1(Eigen::MatrixXd& mat);
  void MultiplyBlocks(const Eigen::Block<const Eigen::MatrixXd>& mat, Index i1,
                      Index i2);

  Eigen::MatrixXd getReductionVar();

 private:
  const Eigen::MatrixXd* rightoperator_ = nullptr;
  const Eigen::MatrixXd* leftoperator_ = nullptr;
  const Eigen::VectorXd* vec_ = nullptr;
  std::vector<Eigen::MatrixXd*> cpu_intermediate_input_;
  std::vector<Eigen::MatrixXd> cpu_temporaries_;

  std::vector<Eigen::MatrixXd> reduction_;
  std::vector<Eigen::VectorXd> temp_;

  bool inside_Parallel_region_;
  Index threadID_parent_;

  Index getParentThreadId() const;

  Index getLocalThreadId(Index ParentThreadId) const;

  Index getNumberThreads() const;

#ifdef USE_CUDA
  bool isGPUthread(Index ParentThreadId) const;

  struct GPU_data {

    explicit GPU_data(Index i)
        : Id(i), pipeline(std::make_unique<CudaPipeline>(int(i))) {
      ;
    }

    Index Id;
    std::unique_ptr<CudaPipeline> pipeline;
    std::vector<std::unique_ptr<CudaMatrix>> temp;

    CudaMatrix& Mat(Index i) { return *temp[i]; }
    CudaPipeline& pipe() { return *pipeline; }
    void activateGPU() { checkCuda(cudaSetDevice(pipeline->getDeviceId())); }

    void push_back(const Eigen::MatrixXd& m) {
      temp.push_back(std::make_unique<CudaMatrix>(m, pipeline->get_stream()));
    }
    void push_back(Index rows, Index cols) {
      temp.push_back(
          std::make_unique<CudaMatrix>(rows, cols, pipeline->get_stream()));
    }

    void resize(Index id, Index rows, Index cols) {
      temp[id] =
          std::make_unique<CudaMatrix>(rows, cols, pipeline->get_stream());
    }
  };

  std::vector<GPU_data> gpus_;
  static bool isInVector(Index Id, const std::vector<GPU_data>& vec);
#endif
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_OPENMP_CUDA_H
