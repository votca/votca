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

#define BOOST_TEST_MODULE symmetric_matrix_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/xtp/symmetric_matrix.h>

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(symmetric_matrix_test)

BOOST_AUTO_TEST_CASE(Constructor_test) {
  Index dim = 3;
  Eigen::MatrixXd test = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd trans = test.transpose();
  test += trans;
  Symmetric_Matrix sym = Symmetric_Matrix(test);

  Eigen::MatrixXcd back = Eigen::MatrixXcd::Zero(dim, dim);

  sym.AddtoEigenMatrix(back);

  bool check_matrices = back.isApprox(test, 0.000001);
  if (!check_matrices) {
    std::cout << sym << std::endl;
    std::cout << test << std::endl;
    std::cout << back << std::endl;
  }
  BOOST_CHECK_EQUAL(check_matrices, 1);
}

BOOST_AUTO_TEST_CASE(Add_test) {

  Index dim = 3;
  Eigen::MatrixXd test = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd trans = test.transpose();
  test += trans;
  Symmetric_Matrix sym = Symmetric_Matrix(test);

  Eigen::MatrixXcd rand = Eigen::MatrixXcd::Random(dim, dim);

  Eigen::MatrixXcd result = rand + 2 * test;

  sym.AddtoEigenMatrix(rand, 2.0);

  bool check_matrices2 = rand.isApprox(result, 0.000001);
  if (!check_matrices2) {
    std::cout << test << std::endl;
    std::cout << sym << std::endl;
    std::cout << rand << std::endl;
    std::cout << result << std::endl;
  }
  BOOST_CHECK_EQUAL(check_matrices2, 1);
}

BOOST_AUTO_TEST_CASE(AddUpper_test) {

  Index dim = 3;
  Eigen::MatrixXd test = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd trans = test.transpose();
  test += trans;
  Symmetric_Matrix sym = Symmetric_Matrix(test);

  Eigen::MatrixXd rand = Eigen::MatrixXd::Random(dim, dim);

  Eigen::MatrixXd result = rand + 2 * test;
  Eigen::SelfAdjointView<Eigen::MatrixXd, Eigen::Upper> m =
      rand.selfadjointView<Eigen::Upper>();
  sym.AddtoEigenUpperMatrix(m, 2.0);

  result.triangularView<Eigen::StrictlyLower>() =
      Eigen::MatrixXd::Zero(dim, dim);
  rand.triangularView<Eigen::StrictlyLower>() = Eigen::MatrixXd::Zero(dim, dim);

  bool check_matrices2 = rand.isApprox(result, 0.000001);
  if (!check_matrices2) {
    std::cout << "ref" << std::endl;
    std::cout << rand << std::endl;
    std::cout << "sym matrix full" << std::endl;
    std::cout << result << std::endl;
  }
  BOOST_CHECK_EQUAL(check_matrices2, 1);
}

BOOST_AUTO_TEST_CASE(FullMatrix_test) {

  Index dim = 3;
  Eigen::MatrixXd test = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd trans = test.transpose();
  test += trans;
  Symmetric_Matrix sym = Symmetric_Matrix(test);
  Eigen::MatrixXd result = sym.FullMatrix();
  bool check_matrices = test.isApprox(result, 0.000001);

  BOOST_CHECK_EQUAL(check_matrices, 1);
}

BOOST_AUTO_TEST_CASE(UpperMatrix_test) {

  Index dim = 3;
  Eigen::MatrixXd test = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd trans = test.transpose();
  test += trans;

  Symmetric_Matrix sym = Symmetric_Matrix(test);
  Eigen::MatrixXd result = sym.UpperMatrix();

  result.triangularView<Eigen::StrictlyLower>() =
      Eigen::MatrixXd::Zero(dim, dim);
  test.triangularView<Eigen::StrictlyLower>() = Eigen::MatrixXd::Zero(dim, dim);

  bool check_matrices = test.isApprox(result, 0.000001);

  BOOST_CHECK_EQUAL(check_matrices, 1);
}

BOOST_AUTO_TEST_CASE(TraceofProd_test) {

  Index dim = 3;
  Eigen::MatrixXd test = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd trans = test.transpose();
  test += trans;
  Symmetric_Matrix sym1 = Symmetric_Matrix(test);

  Eigen::MatrixXd test2 = Eigen::MatrixXd::Random(dim, dim);
  Eigen::MatrixXd trans2 = test2.transpose();
  test2 += trans2;
  Symmetric_Matrix sym2 = Symmetric_Matrix(test2);

  double ref = test2.cwiseProduct(test).sum();
  double result = sym1.TraceofProd(sym2);

  bool check = std::abs(ref - result) < 1e-7;
  if (!check) {
    std::cout << "Ref: " << ref << " Sym: " << result << std::endl;
  }

  BOOST_CHECK_EQUAL(check, 1);
}

BOOST_AUTO_TEST_SUITE_END()
