/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>
#include <exception>
#include <iostream>
#include <votca/tools/table.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(table_test)

BOOST_AUTO_TEST_CASE(create_test) { Table tb; }

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

BOOST_AUTO_TEST_CASE(xy_test) {
  Table tb;
  for (double x = 0; x < 10; ++x) {
    double y = 2 * x;
    tb.push_back(x, y);
  }

  auto x_v = tb.x();
  auto y_v = tb.y();
  for (int i = 0; i < 10; ++i) {
    int x = i;
    int y = 2 * x;
    BOOST_CHECK_EQUAL(boost::lexical_cast<int>(x_v(i)),
                      boost::lexical_cast<int>(tb.x(i)));
    BOOST_CHECK_EQUAL(boost::lexical_cast<int>(y_v(i)),
                      boost::lexical_cast<int>(tb.y(i)));
    BOOST_CHECK_EQUAL(x, boost::lexical_cast<int>(tb.x(i)));
    BOOST_CHECK_EQUAL(y, boost::lexical_cast<int>(tb.y(i)));
    BOOST_CHECK_EQUAL(x, boost::lexical_cast<int>(x_v(i)));
    BOOST_CHECK_EQUAL(y, boost::lexical_cast<int>(y_v(i)));
  }
}

BOOST_AUTO_TEST_CASE(getMinMax_test) {
  Table tb;
  for (double x = 0; x < 10; ++x) {
    double y = 2 * x;
    tb.push_back(x, y);
  }
  
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(tb.getMinX()), 0);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(tb.getMaxX()), 9);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(tb.getMinY()), 0);
  BOOST_CHECK_EQUAL(boost::lexical_cast<int>(tb.getMaxY()), 18);
}

BOOST_AUTO_TEST_SUITE_END()
