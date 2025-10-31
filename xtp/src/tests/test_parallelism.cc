/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE parallelism_test

// Third party includes
#include <boost/test/unit_test.hpp>

#include "votca/xtp/eigen.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(parallelism_test)
#if defined(_OPENMP)
BOOST_AUTO_TEST_CASE(openmp) {

  // Check if number is set correctly
  int inThreads = 2;
  OPENMP::setMaxThreads(inThreads);
  int outThreads = static_cast<int>(OPENMP::getMaxThreads());
  BOOST_CHECK(inThreads == outThreads);

  // Check if the number of threads is actually generated
  int sumThreads = 0;
#pragma omp parallel reduction(+ : sumThreads)
  {
    sumThreads += 1;
  }
  BOOST_CHECK(inThreads == sumThreads);
}
#endif
BOOST_AUTO_TEST_SUITE_END()
