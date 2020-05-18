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

#define BOOST_TEST_MODULE histogramnew_test
#include "../../include/votca/tools/histogramnew.h"
#include "../../include/votca/tools/table.h"
#include <boost/test/unit_test.hpp>
#include <exception>
#include <iostream>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(histogramnew_test)

BOOST_AUTO_TEST_CASE(create_test) { HistogramNew hn; }

BOOST_AUTO_TEST_CASE(init_test) {
  HistogramNew hn;
  double min_v = 1.2;
  double max_v = 102.0;
  hn.Initialize(min_v, max_v, 10);
}

BOOST_AUTO_TEST_CASE(step_test) {
  HistogramNew hn;
  double min_v = 1;
  double max_v = 9;
  hn.Initialize(min_v, max_v, 8);
  double step = hn.getStep();
  votca::Index value = static_cast<votca::Index>(step);
  BOOST_CHECK_EQUAL(value * 10, 10);
}

BOOST_AUTO_TEST_CASE(nbins_test) {
  HistogramNew hn;
  double min_v = 1.0;
  double max_v = 9.0;
  hn.Initialize(min_v, max_v, 8);
  votca::Index bins = hn.getNBins();
  BOOST_CHECK_EQUAL(bins * 10, 80);
}

BOOST_AUTO_TEST_CASE(Process_test) {
  HistogramNew hn;
  double min_v = 0.0;
  double max_v = 10.0;
  hn.Initialize(min_v, max_v, 11);
  vector<double> data;
  for (double x = 0; x < 10; ++x) {
    data.push_back(x);
  }
  hn.ProcessRange(data.begin(), data.end());
  hn.Process(4.5);
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(hn.getStep() * 10), 10);
  auto dat = hn.data();
  // Range -0.5 - 0.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(0)), 1);
  // Range 0.5 - 1.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(1)), 1);
  // Range 1.5 - 2.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(2)), 1);
  // Range 2.5 - 3.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(3)), 1);
  // Range 3.5 - 4.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(4)), 1);
  // Range 4.5 - 5.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(5)), 2);
  // Range 5.5 - 6.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(6)), 1);
  // Range 6.5 - 7.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(7)), 1);
  // Range 7.5 - 8.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(8)), 1);
  // Range 8.5 - 9.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(9)), 1);
  // Range 9.5 - 10.5
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(dat.y(10)), 0);
}

BOOST_AUTO_TEST_CASE(minmax_test) {
  HistogramNew hn;
  double min_v = 0.0;
  double max_v = 10.0;
  hn.Initialize(min_v, max_v, 10);
  vector<double> data;
  for (double x = 0; x < 9; ++x) {
    data.push_back(x);
  }
  hn.ProcessRange(data.begin(), data.end());
  hn.Process(4.5);
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(hn.getMinBinVal()), 0);
  BOOST_CHECK_EQUAL(static_cast<votca::Index>(hn.getMaxBinVal()), 2);
}

BOOST_AUTO_TEST_SUITE_END()
