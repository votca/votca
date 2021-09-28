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
#define BOOST_TEST_MODULE CudaMatrix_test

// Standard includes
#include <string>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/cudapipeline.h"
#include "votca/xtp/eigen.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(CudaMatrix_test)

BOOST_AUTO_TEST_CASE(create_cudamatrix) {
  // Call the class to handle GPU resources
  CudaPipeline cp(0);

  // Call matrix multiplication GPU
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(5, 7);

  // Copy matrix back and for to the GPU
  CudaMatrix cumatrix{B, cp.get_stream()};
  Eigen::MatrixXd tmp = cumatrix;

  // Expected results
  BOOST_TEST(B.isApprox(tmp));
}

BOOST_AUTO_TEST_CASE(create_cudamatrixblock) {
  // Call the class to handle GPU resources
  CudaPipeline cp(0);

  // Call matrix multiplication GPU
  Eigen::MatrixXd X = Eigen::MatrixXd::Random(10, 8);

  // Copy matrix back and for to the GPU
  CudaMatrix cumatrix{X.block(2, 4, 3, 1), cp.get_stream()};
  Eigen::MatrixXd tmp = cumatrix;

  // Expected results
  BOOST_TEST(X.block(2, 4, 3, 1).isApprox(tmp));
}

BOOST_AUTO_TEST_CASE(setZero) {
  // Call the class to handle GPU resources
  CudaPipeline cp(0);

  // Call matrix multiplication GPU
  Eigen::MatrixXd X = Eigen::MatrixXd::Random(10, 8);

  // Copy matrix back and for to the GPU
  CudaMatrix cumatrix{X, cp.get_stream()};
  cumatrix.setZero();
  Eigen::MatrixXd tmp = cumatrix;

  // Expected results
  BOOST_TEST(tmp.isApproxToConstant(0, 1e-14));
}

BOOST_AUTO_TEST_SUITE_END()
