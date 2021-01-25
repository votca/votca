/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE CudaPipeline_test

// Standard includes
#include <string>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/cudapipeline.h"
#include "votca/xtp/eigen.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(CudaPipeline_test)

BOOST_AUTO_TEST_CASE(matmul) {
  // Call the class to handle GPU resources
  CudaPipeline cuda_pip(0);

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(6, 10);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(10, 6);

  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6, 6);

  CudaMatrix Ag{A, cuda_pip.get_stream()};
  CudaMatrix Bg{B, cuda_pip.get_stream()};
  CudaMatrix Cg{C, cuda_pip.get_stream()};

  cuda_pip.gemm(Ag, Bg, Cg);

  Eigen::MatrixXd GPU_result = Cg;
  Eigen::MatrixXd CPU_result = A * B;
  bool check = CPU_result.isApprox(GPU_result, 1e-9);
  BOOST_CHECK_EQUAL(check, true);
  if (!check) {
    std::cout << "CPU\n" << CPU_result << std::endl;
    std::cout << "GPU\n" << GPU_result << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(matmul_add) {
  // Call the class to handle GPU resources
  CudaPipeline cuda_pip(0);

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(6, 10);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(10, 6);

  Eigen::MatrixXd C = Eigen::MatrixXd::Random(6, 6);

  CudaMatrix Ag{A, cuda_pip.get_stream()};
  CudaMatrix Bg{B, cuda_pip.get_stream()};
  CudaMatrix Cg{C, cuda_pip.get_stream()};

  cuda_pip.gemm(Ag, Bg, Cg, 1.0);

  Eigen::MatrixXd GPU_result = Cg;
  Eigen::MatrixXd CPU_result = A * B + C;
  bool check = CPU_result.isApprox(GPU_result, 1e-9);
  BOOST_CHECK_EQUAL(check, true);
  if (!check) {
    std::cout << "CPU\n" << CPU_result << std::endl;
    std::cout << "GPU\n" << GPU_result << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(matmul_ABt) {
  // Call the class to handle GPU resources
  CudaPipeline cuda_pip(0);

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(6, 10);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(6, 10);

  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6, 6);

  CudaMatrix Ag{A, cuda_pip.get_stream()};
  CudaMatrix Bg{B, cuda_pip.get_stream()};
  CudaMatrix Cg{C, cuda_pip.get_stream()};

  cuda_pip.gemm(Ag, Bg.transpose(), Cg);

  Eigen::MatrixXd GPU_result = Cg;
  Eigen::MatrixXd CPU_result = A * B.transpose();
  bool check = CPU_result.isApprox(GPU_result, 1e-9);
  BOOST_CHECK_EQUAL(check, true);
  if (!check) {
    std::cout << "CPU\n" << CPU_result << std::endl;
    std::cout << "GPU\n" << GPU_result << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(matmul_AtB) {
  // Call the class to handle GPU resources
  CudaPipeline cuda_pip(0);

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 6);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(10, 6);

  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6, 6);

  CudaMatrix Ag{A, cuda_pip.get_stream()};
  CudaMatrix Bg{B, cuda_pip.get_stream()};
  CudaMatrix Cg{C, cuda_pip.get_stream()};

  cuda_pip.gemm(Ag.transpose(), Bg, Cg);

  Eigen::MatrixXd GPU_result = Cg;
  Eigen::MatrixXd CPU_result = A.transpose() * B;
  bool check = CPU_result.isApprox(GPU_result, 1e-9);
  BOOST_CHECK_EQUAL(check, true);
  if (!check) {
    std::cout << "CPU\n" << CPU_result << std::endl;
    std::cout << "GPU\n" << GPU_result << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(matmul_AtBt) {
  // Call the class to handle GPU resources
  CudaPipeline cuda_pip(0);

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 6);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(6, 10);

  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6, 6);

  CudaMatrix Ag{A, cuda_pip.get_stream()};
  CudaMatrix Bg{B, cuda_pip.get_stream()};
  CudaMatrix Cg{C, cuda_pip.get_stream()};

  cuda_pip.gemm(Ag.transpose(), Bg.transpose(), Cg);

  Eigen::MatrixXd GPU_result = Cg;
  Eigen::MatrixXd CPU_result = A.transpose() * B.transpose();
  bool check = CPU_result.isApprox(GPU_result, 1e-9);
  BOOST_CHECK_EQUAL(check, true);
  if (!check) {
    std::cout << "CPU\n" << CPU_result << std::endl;
    std::cout << "GPU\n" << GPU_result << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(diag_matrix_mul) {
  // Call the class to handle GPU resources
  CudaPipeline cuda_pip(0);

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 6);
  Eigen::VectorXd b = Eigen::VectorXd::Random(10);

  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(10, 6);

  CudaMatrix Ag{A, cuda_pip.get_stream()};
  CudaMatrix bg{b, cuda_pip.get_stream()};
  CudaMatrix Cg{C, cuda_pip.get_stream()};

  cuda_pip.diag_gemm(Ag, bg, Cg);

  Eigen::MatrixXd GPU_result = Cg;
  Eigen::MatrixXd CPU_result = b.asDiagonal() * A;
  bool check = CPU_result.isApprox(GPU_result, 1e-9);
  BOOST_CHECK_EQUAL(check, true);
  if (!check) {
    std::cout << "CPU\n" << CPU_result << std::endl;
    std::cout << "GPU\n" << GPU_result << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
