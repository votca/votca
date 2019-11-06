/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE cubicspline_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/tools/cubicspline.h>

using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(cubicspline_test)

BOOST_AUTO_TEST_CASE(cubicspline_fit_test) {

  votca::Index size = 80;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd y = Eigen::VectorXd::Zero(size);
  for (votca::Index i = 0; i < size; ++i) {
    x(i) = 0.25 * double(i);
    y(i) = std::sin(x(i));
  }
  CubicSpline cspline;
  cspline.setBCInt(0);
  cspline.GenerateGrid(0.4, 0.6, 0.1);
  cspline.Fit(x, y);

  Eigen::VectorXd rs = Eigen::VectorXd::Zero(10);
  rs << 0.45, 0.47, 0.8, 0.75, 0.6, 0.4, 0.9, 0.55, 0, 0;
  Eigen::VectorXd values_ref = Eigen::VectorXd::Zero(10);
  values_ref << 0.311213, 0.310352, 0.296153, 0.298304, 0.304759, 0.313364,
      0.291851, 0.30691, 0.33058, 0.33058;
  Eigen::VectorXd derivatives_ref = Eigen::VectorXd::Zero(10);
  derivatives_ref << -0.0430277, -0.0430282, -0.0430231, -0.0430267, -0.0430313,
      -0.0430272, -0.0430128, -0.0430308, -0.04306, -0.04306;
  Eigen::VectorXd values = cspline.Calculate(rs);
  Eigen::VectorXd derivatives = cspline.CalculateDerivative(rs);

  bool equal_val = values_ref.isApprox(values, 1e-5);

  if (!equal_val) {
    std::cout << "result value" << std::endl;
    std::cout << values.transpose() << std::endl;
    std::cout << "ref value" << std::endl;
    std::cout << values_ref.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equal_val, true);

  bool equal_derivative = derivatives_ref.isApprox(derivatives, 1e-5);

  if (!equal_derivative) {
    std::cout << "result value" << std::endl;
    std::cout << derivatives.transpose() << std::endl;
    std::cout << "ref value" << std::endl;
    std::cout << derivatives_ref.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equal_derivative, true);
}

BOOST_AUTO_TEST_CASE(cubicspline_interpolate_test) {

  int size = 80;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd y = Eigen::VectorXd::Zero(size);
  for (int i = 0; i < size; ++i) {
    x(i) = 0.25 * i;
    y(i) = std::sin(x(i));
  }
  CubicSpline cspline;
  cspline.setBCInt(0);
  cspline.Interpolate(x, y);

  Eigen::VectorXd rs = Eigen::VectorXd::Zero(10);
  rs << 0.45, 0.47, 0.8, 0.75, 0.6, 0.4, 0.9, 0.55, 0, 0;
  Eigen::VectorXd values_ref = Eigen::VectorXd::Zero(10);
  values_ref << 0.434964, 0.452886, 0.717353, 0.681639, 0.564637, 0.389415,
      0.78332, 0.522684, 0, 0;
  Eigen::VectorXd derivatives_ref = Eigen::VectorXd::Zero(10);
  derivatives_ref << 0.900494, 0.891601, 0.69661, 0.731673, 0.825305, 0.921091,
      0.621663, 0.852451, 0.999978, 0.999978;
  Eigen::VectorXd values = cspline.Calculate(rs);
  Eigen::VectorXd derivatives = cspline.CalculateDerivative(rs);

  bool equal_val = values_ref.isApprox(values, 1e-5);

  if (!equal_val) {
    std::cout << "result value" << std::endl;
    std::cout << values.transpose() << std::endl;
    std::cout << "ref value" << std::endl;
    std::cout << values_ref.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equal_val, true);

  bool equal_derivative = derivatives_ref.isApprox(derivatives, 1e-5);

  if (!equal_derivative) {
    std::cout << "result value" << std::endl;
    std::cout << derivatives.transpose() << std::endl;
    std::cout << "ref value" << std::endl;
    std::cout << derivatives_ref.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equal_derivative, true);
}

BOOST_AUTO_TEST_CASE(cubicspline_matrix_test) {

  CubicSpline cspline;
  cspline.setBCInt(0);
  cspline.GenerateGrid(0.4, 0.6, 0.1);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(1, 6);
  Eigen::MatrixXd Aref = Eigen::MatrixXd::Zero(1, 6);

  Aref << 0.0, -9.0, 10.0, 0.0, -0.03333333, -0.01666667;

  cspline.AddToFitMatrix(A, 0.5, 0, 0, 1.0, 1.0);

  bool equalMatrix = Aref.isApprox(A, 1e-5);
  if (!equalMatrix) {
    std::cout << "result A" << std::endl;
    std::cout << A << std::endl;
    std::cout << "ref A" << std::endl;
    std::cout << Aref << std::endl;
  }
  BOOST_CHECK_EQUAL(equalMatrix, true);
}

BOOST_AUTO_TEST_SUITE_END()
