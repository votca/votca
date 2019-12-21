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

#define BOOST_TEST_MODULE floatingpointcomparison_test
#include "../../include/votca/tools/floatingpointcomparison.h"
#include <boost/test/unit_test.hpp>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(floatingpointcomparison_test)

BOOST_AUTO_TEST_CASE(comparison_test) {
  double var1 = 1.000;
  double var2 = 1.010;
  BOOST_CHECK(isApproximatelyEqual(var1, var2, 0.1));
  BOOST_CHECK(isApproximatelyEqual(var1, var2, 0.01));
  BOOST_CHECK(!isApproximatelyEqual(var1, var2, 0.001));
}

BOOST_AUTO_TEST_SUITE_END()
