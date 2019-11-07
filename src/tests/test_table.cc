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

#define BOOST_TEST_MODULE table_test
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include <exception>
#include <iostream>
#include <votca/tools/table.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(table_test)

BOOST_AUTO_TEST_CASE(create_test) {
  Table tb;
  Table tb2(tb);
}

BOOST_AUTO_TEST_CASE(size_test) {
  Table tb;
  BOOST_CHECK_EQUAL(tb.size(), 0);
}

BOOST_AUTO_TEST_CASE(pushback_test) {
  Table tb;
  for (double x = 0; x < 10; ++x) {
    double y = 2 * x;
    tb.push_back(x, y);
  }
  BOOST_CHECK_EQUAL(tb.size(), 10);
}

BOOST_AUTO_TEST_CASE(resize_test) {
  Table tb;
  tb.resize(0);
  tb.resize(10);

  bool error_thrown = false;
  try {
    tb.resize(-5);
  } catch (...) {
    error_thrown = true;
  }
  BOOST_CHECK(error_thrown);
}

BOOST_AUTO_TEST_CASE(xy_test) {
  Table tb;
  for (double x = 0; x < 10; ++x) {
    double y = 2 * x;
    tb.push_back(x, y);
  }

  auto x_v = tb.x();
  auto y_v = tb.y();
  for (votca::Index i = 0; i < 10; ++i) {
    votca::Index x = i;
    votca::Index y = 2 * x;
    BOOST_CHECK_EQUAL(static_cast<votca::Index>(x_v(i)),
                      static_cast<votca::Index>(tb.x(i)));
    BOOST_CHECK_EQUAL(static_cast<votca::Index>(y_v(i)),
                      static_cast<votca::Index>(tb.y(i)));
    BOOST_CHECK_EQUAL(x, static_cast<votca::Index>(tb.x(i)));
    BOOST_CHECK_EQUAL(y, static_cast<votca::Index>(tb.y(i)));
    BOOST_CHECK_EQUAL(x, static_cast<votca::Index>(x_v(i)));
    BOOST_CHECK_EQUAL(y, static_cast<votca::Index>(y_v(i)));
  }
}

BOOST_AUTO_TEST_CASE(getMinMax_test) {
  Table tb;
  for (double x = 0; x < 10; ++x) {
    double y = 2 * x;
    tb.push_back(x, y);
  }

  BOOST_CHECK_EQUAL(static_cast<votca::Index>(tb.getMinX()), 0);
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(tb.getMaxX()), 9);
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(tb.getMinY()), 0);
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(tb.getMaxY()), 18);
}

BOOST_AUTO_TEST_CASE(generate_grid_spacing_test) {
  Table tb;
  double min_v = 1.2;
  double max_v = 2.0;

  tb.GenerateGridSpacing(min_v, max_v, 0.2);

  BOOST_CHECK_EQUAL(tb.size(), 5);
  BOOST_CHECK_CLOSE(tb.getMinX(), 1.2, 1e-5);
  BOOST_CHECK_CLOSE(tb.getMaxX(), 2.0, 1e-5);
}

BOOST_AUTO_TEST_CASE(smoothing_test) {
  Table tb;
  double min_v = 1.2;
  double max_v = 2.0;

  tb.GenerateGridSpacing(min_v, max_v, 0.1);

  BOOST_CHECK_EQUAL(tb.size(), 9);
  tb.y() = tb.x().array().sinh();
  tb.Smooth(2);
  Eigen::VectorXd refy = Eigen::VectorXd::Zero(9);
  refy << 1.50946, 1.70621, 1.91527, 2.14055, 2.38962, 2.65963, 2.95959,
      3.28268, 3.62686;

  bool equal = tb.y().isApprox(refy, 1e-5);

  if (!equal) {
    std::cout << "result value" << std::endl;
    std::cout << tb.y().transpose() << std::endl;
    std::cout << "ref value" << std::endl;
    std::cout << refy.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equal, true);
}

BOOST_AUTO_TEST_SUITE_END()
