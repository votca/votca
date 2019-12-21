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

#define BOOST_TEST_MODULE linspline_test
#include "../../include/votca/tools/linspline.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(linspline_test)

BOOST_AUTO_TEST_CASE(fit_test) {

  int size = 80;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd y = Eigen::VectorXd::Zero(size);
  for (int i = 0; i < size; ++i) {
    x(i) = 0.25 * i;
    y(i) = std::sin(x(i));
  }
  LinSpline cspline;
  cspline.setBCInt(0);
  cspline.GenerateGrid(0.4, 0.6, 0.1);
  cspline.Fit(x, y);

  Eigen::VectorXd rs = Eigen::VectorXd::Zero(10);
  rs << 0.45, 0.47, 0.8, 0.75, 0.6, 0.4, 0.9, 0.55, 0, 0;
  Eigen::VectorXd values_ref = Eigen::VectorXd::Zero(10);
  values_ref << 0.218354, 0.225766, 0.230181, 0.231299, 0.23465, 0.199824,
      0.227947, 0.235767, 0.0515847, 0.0515847;
  Eigen::VectorXd derivatives_ref = Eigen::VectorXd::Zero(10);
  derivatives_ref << 0.370599, 0.370599, -0.0223429, -0.0223429, -0.0223429,
      0.370599, -0.0223429, -0.0223429, 0.370599, 0.370599;
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

BOOST_AUTO_TEST_CASE(interpolate_test) {

  int size = 80;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd y = Eigen::VectorXd::Zero(size);
  for (int i = 0; i < size; ++i) {
    x(i) = 0.25 * i;
    y(i) = std::sin(x(i));
  }
  LinSpline cspline;
  cspline.setBCInt(0);
  cspline.Interpolate(x, y);

  Eigen::VectorXd rs = Eigen::VectorXd::Zero(10);
  rs << 0.45, 0.47, 0.8, 0.75, 0.6, 0.4, 0.9, 0.55, 0, 0;
  Eigen::VectorXd values_ref = Eigen::VectorXd::Zero(10);
  values_ref << 0.433021, 0.451583, 0.713605, 0.681639, 0.560311, 0.386617,
      0.777538, 0.519868, 0, 0;
  Eigen::VectorXd derivatives_ref = Eigen::VectorXd::Zero(10);
  derivatives_ref << 0.928086, 0.928086, 0.639329, 0.639329, 0.808853, 0.928086,
      0.639329, 0.808853, 0.989616, 0.989616;
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

BOOST_AUTO_TEST_SUITE_END()
