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

#define BOOST_TEST_MODULE edge_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <set>
#include <unordered_map>

#include <votca/tools/edge.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(edge_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Edge ed(2, 3); }

BOOST_AUTO_TEST_CASE(ouput_test) {
  Edge ed(2, 3);
  cout << ed << endl;
}

BOOST_AUTO_TEST_CASE(equivalence_test) {
  Edge ed(2, 3);
  Edge ed2(2, 4);
  Edge ed3(3, 2);
  BOOST_CHECK_EQUAL(ed, ed);
  BOOST_CHECK_EQUAL(ed, ed3);
  BOOST_CHECK_EQUAL((ed == ed2),false);
}

BOOST_AUTO_TEST_CASE(nequivalence_test) {
  Edge ed(2, 3);
  Edge ed2(2, 4);
  BOOST_CHECK_EQUAL(ed != ed2, true);
  BOOST_CHECK_EQUAL((ed != ed), false);
}

BOOST_AUTO_TEST_CASE(getter_test) {
  Edge ed(2, 3);
  BOOST_CHECK_EQUAL(ed.getV1(), 2);
  BOOST_CHECK_EQUAL(ed.getV2(), 3);
  Edge ed2(3,2);
  BOOST_CHECK_EQUAL(ed.getV1(), 2);
  BOOST_CHECK_EQUAL(ed.getV2(), 3);
}

BOOST_AUTO_TEST_CASE(less_test) {
  Edge ed1(1, 2);
  Edge ed2(2, 1);
  BOOST_CHECK_EQUAL(!(ed1 < ed2), true);
  Edge ed3(3, 0);
  BOOST_CHECK_EQUAL(ed3 < ed1, true);
  Edge ed4(1, 5);
  BOOST_CHECK_EQUAL(ed1 < ed4, true);
  Edge ed5(2, 2);
  BOOST_CHECK_EQUAL(ed4 < ed5, true);
}

BOOST_AUTO_TEST_CASE(greater_test) {
  Edge ed1(1, 2);
  Edge ed2(2, 1);
  BOOST_CHECK_EQUAL(!(ed1 > ed2), true);
  Edge ed3(3, 0);
  BOOST_CHECK_EQUAL(!(ed3 > ed1), true);
}

BOOST_AUTO_TEST_CASE(hash_key_test) {
  set<Edge> e_set;
  Edge      ed(23, 43);
  e_set.insert(ed);
  unordered_map<int, Edge> e_map;
  e_map[2] = ed;
}

BOOST_AUTO_TEST_SUITE_END()
