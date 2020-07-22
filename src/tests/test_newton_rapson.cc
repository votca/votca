/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE newton_rapson

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/newton_rapson.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(newton_rapson)
class Func {
 public:
  std::pair<double, double> operator()(double x) const {
    std::pair<double, double> value;
    value.first = x * x - 612;
    value.second = 2 * x;
    return value;
  }
};

BOOST_AUTO_TEST_CASE(newton_rapson) {

  Func f;
  Index iterations = 50;
  double tolerance = 1e-9;
  NewtonRapson<Func> n(iterations, tolerance);
  double root = n.FindRoot(f, 10.0);
  BOOST_CHECK_EQUAL(n.getInfo(), NewtonRapson<Func>::Errors::success);
  // see https://en.wikipedia.org/wiki/Newton%27s_method
  BOOST_CHECK_CLOSE(root, 24.738633753, 1e-7);
}

BOOST_AUTO_TEST_SUITE_END()
