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

#define BOOST_TEST_MODULE edge_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <set>
#include <unordered_map>

#include <votca/tools/reducededge.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(edge_test)

BOOST_AUTO_TEST_CASE(constructors_test) { ReducedEdge ed(2, 3); }

BOOST_AUTO_TEST_CASE(ouput_test) {
  ReducedEdge ed(2, 3);
  cout << ed << endl;
}

BOOST_AUTO_TEST_CASE(equivalence_test) {
  ReducedEdge ed(2, 3);
  ReducedEdge ed2(2, 4);
  ReducedEdge ed3(3, 2);
  BOOST_CHECK_EQUAL(ed, ed);
  BOOST_CHECK_EQUAL(ed, ed3);
  BOOST_CHECK_EQUAL((ed == ed2), false);

  ReducedEdge ed4(vector<votca::Index>{2, 3, 4});
  ReducedEdge ed5(vector<votca::Index>{4, 3, 2});
  BOOST_CHECK_EQUAL(ed4, ed5);
  BOOST_CHECK_EQUAL((ed4 == ed), false);

  ReducedEdge ed6(vector<votca::Index>{2, 6, 3, 1});
  ReducedEdge ed7(vector<votca::Index>{1, 3, 6, 2});

  auto chain = ed6.getChain();
  BOOST_CHECK_EQUAL(ed6, ed7);

  // Both are loops they get rearranged internally so the smallest vertex
  // appears as EndPoint1 and EndPoint2
  ReducedEdge ed8(vector<votca::Index>{1, 2, 4, 1});
  ReducedEdge ed9(vector<votca::Index>{2, 4, 1, 2});
  BOOST_CHECK_EQUAL(ed8, ed9);
}

BOOST_AUTO_TEST_CASE(nequivalence_test) {
  ReducedEdge ed(2, 3);
  ReducedEdge ed2(2, 4);
  BOOST_CHECK_EQUAL(ed != ed2, true);
  BOOST_CHECK_EQUAL((ed != ed), false);
}

BOOST_AUTO_TEST_CASE(getter_test) {
  ReducedEdge ed(2, 3);
  BOOST_CHECK_EQUAL(ed.getEndPoint1(), 2);
  BOOST_CHECK_EQUAL(ed.getEndPoint2(), 3);
  ReducedEdge ed2(3, 2);
  BOOST_CHECK_EQUAL(ed.getEndPoint1(), 2);
  BOOST_CHECK_EQUAL(ed.getEndPoint2(), 3);

  ReducedEdge ed3(vector<votca::Index>{2, 6, 3, 1});
  BOOST_CHECK_EQUAL(ed3.getEndPoint1(), 1);
  BOOST_CHECK_EQUAL(ed3.getEndPoint2(), 2);
  ReducedEdge ed4(vector<votca::Index>{1, 3, 6, 2});
  BOOST_CHECK_EQUAL(ed4.getEndPoint1(), 1);
  BOOST_CHECK_EQUAL(ed4.getEndPoint2(), 2);

  ReducedEdge ed5(vector<votca::Index>{1, 2, 4, 1});
  BOOST_CHECK_EQUAL(ed5.getEndPoint1(), 1);
  BOOST_CHECK_EQUAL(ed5.getEndPoint2(), 1);
  ReducedEdge ed6(vector<votca::Index>{2, 4, 1, 2});
  BOOST_CHECK_EQUAL(ed6.getEndPoint1(), 1);
  BOOST_CHECK_EQUAL(ed6.getEndPoint2(), 1);
}

BOOST_AUTO_TEST_CASE(getchain_test) {
  vector<votca::Index> vertices{4, 1, 6, 7, 0};
  ReducedEdge ed(vertices);
  auto chain = ed.getChain();
  BOOST_CHECK_EQUAL(chain.at(0), 0);
  BOOST_CHECK_EQUAL(chain.at(1), 7);
  BOOST_CHECK_EQUAL(chain.at(2), 6);
  BOOST_CHECK_EQUAL(chain.at(3), 1);
  BOOST_CHECK_EQUAL(chain.at(4), 4);
}

BOOST_AUTO_TEST_CASE(vertexwithinchain_test) {
  vector<votca::Index> vertices{1, 3, 4, 5};
  ReducedEdge ed(vertices);
  BOOST_CHECK(ed.vertexExistInChain(1));
  BOOST_CHECK(ed.vertexExistInChain(2) == false);
  BOOST_CHECK(ed.vertexExistInChain(3));
  BOOST_CHECK(ed.vertexExistInChain(4));
  BOOST_CHECK(ed.vertexExistInChain(5));
}

BOOST_AUTO_TEST_CASE(less_test) {
  ReducedEdge ed1(1, 2);
  ReducedEdge ed2(2, 1);
  BOOST_CHECK_EQUAL((ed1 < ed2), false);
  ReducedEdge ed3(3, 0);
  BOOST_CHECK_EQUAL(ed3 < ed1, true);
  ReducedEdge ed4(1, 5);
  BOOST_CHECK_EQUAL(ed1 < ed4, true);
  ReducedEdge ed5(2, 2);
  BOOST_CHECK_EQUAL(ed4 < ed5, true);

  ReducedEdge ed6(vector<votca::Index>{1, 2, 3});
  BOOST_CHECK_EQUAL(ed1 < ed6, true);

  ReducedEdge ed7(vector<votca::Index>{1, 2, 3, 4, 1});
  ReducedEdge ed8(vector<votca::Index>{1, 4, 3, 2, 1});
  BOOST_CHECK_EQUAL((ed7 < ed8), false);  // Equal

  ReducedEdge ed9(vector<votca::Index>{1, 4, 5, 2, 1});
  BOOST_CHECK_EQUAL((ed7 < ed9), true);

  ReducedEdge ed10(vector<votca::Index>{1, 4, 5, 0, 1});
  BOOST_CHECK_EQUAL((ed10 < ed9), true);

  ReducedEdge ed11(vector<votca::Index>{1, 0, 5, 2, 1});
  BOOST_CHECK_EQUAL((ed11 < ed9), true);
}

BOOST_AUTO_TEST_CASE(greater_test) {
  ReducedEdge ed1(1, 2);
  ReducedEdge ed2(2, 1);
  BOOST_CHECK_EQUAL(!(ed1 > ed2), true);
  ReducedEdge ed3(3, 0);
  BOOST_CHECK_EQUAL(!(ed3 > ed1), true);
}

BOOST_AUTO_TEST_CASE(hash_key_test) {
  set<ReducedEdge> e_set;
  ReducedEdge ed(23, 43);
  e_set.insert(ed);
  unordered_map<votca::Index, ReducedEdge> e_map;
  e_map[2] = ed;
}

BOOST_AUTO_TEST_SUITE_END()
