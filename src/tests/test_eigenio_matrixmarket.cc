/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

#define BOOST_TEST_MODULE eigenio_matrixmarket

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/eigenio_matrixmarket.h"
#include "votca/tools/votca_tools_config.h"

using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(eigenio_matrixmarket)

BOOST_AUTO_TEST_CASE(readvector_test) {

  Eigen::VectorXd ref = Eigen::VectorXd::Zero(4);
  ref << 1.0, 2.0, 3.0, 4.0;
  std::cout << std::string(TOOLS_CMAKE_BINARY_DIR) << std::endl;
  Eigen::VectorXd readin = EigenIO_MatrixMarket::ReadVector(
      std::string(TOOLS_CMAKE_BINARY_DIR) +
      "/tools/src/tests/DataFiles/eigenio_matrixmarket/eigen_vector.mm");

  bool check = ref.isApprox(readin, 1e-5);

  if (!check) {
    std::cout << "ref" << std::endl;
    std::cout << ref.transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << readin.transpose() << std::endl;
  }
  BOOST_CHECK(check);
}

BOOST_AUTO_TEST_CASE(writevector_test) {

  Eigen::VectorXd test = Eigen::VectorXd::Random(5);
  EigenIO_MatrixMarket::WriteVector("VectorRandom.mm", test);
  Eigen::VectorXd readin = EigenIO_MatrixMarket::ReadVector("VectorRandom.mm");
  bool check = test.isApprox(readin, 1e-5);

  if (!check) {
    std::cout << "ref" << std::endl;
    std::cout << test.transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << readin.transpose() << std::endl;
  }
  BOOST_CHECK(check);
}

BOOST_AUTO_TEST_CASE(readmatrix_test) {

  Eigen::MatrixXd ref = Eigen::MatrixXd::Zero(4, 3);
  ref << 1.0, 5.0, 9.0, 2.0, 6.0, 10.0, 3.0, 7.0, 11.0, 4.0, 8.0, 12.0;

  Eigen::MatrixXd readin = EigenIO_MatrixMarket::ReadMatrix(
      std::string(TOOLS_CMAKE_BINARY_DIR) +
      "/tools/src/tests/DataFiles/eigenio_matrixmarket/eigen_matrix.mm");

  bool check = ref.isApprox(readin, 1e-5);

  if (!check) {
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << readin << std::endl;
  }
  BOOST_CHECK(check);
}

BOOST_AUTO_TEST_CASE(writematrix_test) {

  Eigen::MatrixXd test = Eigen::MatrixXd::Random(4, 3);
  EigenIO_MatrixMarket::WriteMatrix("MatrixRandom.mm", test);
  Eigen::MatrixXd readin = EigenIO_MatrixMarket::ReadMatrix("MatrixRandom.mm");
  bool check = test.isApprox(readin, 1e-5);

  if (!check) {
    std::cout << "ref" << std::endl;
    std::cout << test << std::endl;
    std::cout << "result" << std::endl;
    std::cout << readin << std::endl;
  }
  BOOST_CHECK(check);
}

BOOST_AUTO_TEST_SUITE_END()
