/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#include <boost/test/unit_test.hpp>
#include <string>
#include <votca/xtp/cudapipeline.h>
#include <votca/xtp/eigen.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(CudaMatrix_test)

BOOST_AUTO_TEST_CASE(create_cudamatrix) {
  // Call the class to handle GPU resources
  CudaPipeline cp;

  // Call matrix multiplication GPU
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 2);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2);
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(3, 2);

  // Define matrices
  A << 1., 2., 3., 4.;
  B << 5., 6., 7., 8., 9., 10.;
  X << 23., 34., 31., 46., 39., 58.;

  // Copy matrix back and for to the GPU
  CudaMatrix cumatrix{B, cp.get_stream()};
  Eigen::MatrixXd tmp = cumatrix;
  Eigen::MatrixXd result = tmp * A;

  // Expected results
  BOOST_TEST(X.isApprox(result));
}
BOOST_AUTO_TEST_SUITE_END()
