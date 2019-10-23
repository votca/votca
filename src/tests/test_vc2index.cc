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

#define BOOST_TEST_MODULE vc2index_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/xtp/vc2index.h>

using namespace votca::xtp;
BOOST_AUTO_TEST_SUITE(vc2index_test)

BOOST_AUTO_TEST_CASE(index_test) {
  int vmin = 0;
  int cmin = 10;
  int ctotal = 10;
  int vtotal = 9;

  vc2index vc = vc2index(vmin, cmin, ctotal);

  int v = 3;
  int c = 12;

  BOOST_CHECK_EQUAL(vc.I(v, c), 32);

  BOOST_CHECK_EQUAL(vc.c(vc.I(v, c)), c);
  BOOST_CHECK_EQUAL(vc.v(vc.I(v, c)), v);

  std::vector<int> indexv;
  std::vector<int> indexc;

  for (int v2 = 0; v2 < vtotal; v2++) {
    for (int c2 = 0; c2 < ctotal; c2++) {
      indexv.push_back(vmin + v2);
      indexc.push_back(cmin + c2);
    }
  }

  for (unsigned j = 0; j < indexv.size(); j++) {
    BOOST_CHECK_EQUAL(indexv[j], vc.v(j));
  }

  for (unsigned j = 0; j < indexc.size(); j++) {
    BOOST_CHECK_EQUAL(indexc[j], vc.c(j));
  }
}

BOOST_AUTO_TEST_SUITE_END()
