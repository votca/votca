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

#include <numeric>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE ndimvector_test

// Standard includes
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/NDimVector.h"

using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(ndimvector_test)

BOOST_AUTO_TEST_CASE(constructor) {

  NDimVector<int, 3> a(1, 3, 4);
  BOOST_CHECK_EQUAL(a.dimension(0), 1);
  BOOST_CHECK_EQUAL(a.dimension(1), 3);
  BOOST_CHECK_EQUAL(a.dimension(2), 4);
  BOOST_CHECK_EQUAL(a.size(), 12);

  a(0, 1, 2) = 2;
  a(0, 1, 3) = 4;
  int sum = std::accumulate(a.begin(), a.end(), 0);
  BOOST_CHECK_EQUAL(sum, 6);

  BOOST_CHECK_EQUAL(a(0, 1, 3), 4);

  NDimVector<std::vector<int>, 3> b(2, 3, 4);
  b(0, 2, 2).push_back(1);
}

BOOST_AUTO_TEST_SUITE_END()
