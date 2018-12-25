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

#define BOOST_TEST_MODULE edgecontainer_test
#include <boost/test/unit_test.hpp>
#include <exception>
#include <iostream>
#include <vector>
#include <votca/tools/edgecontainer.h>
#include <votca/tools/edge.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(edgecontainer_test)

BOOST_AUTO_TEST_CASE(create_test) {
  EdgeContainer edCo;
  Edge ed(1, 2);
  EdgeContainer edCo2(ed);
  Edge ed2(3, 4);
  vector<Edge> eds{ed, ed2};
  EdgeContainer edCo3(eds);
}

BOOST_AUTO_TEST_CASE(getdegree_test){
  Edge ed(1, 2);
  Edge ed2(2, 3);
  vector<Edge> eds{ed, ed2};
  EdgeContainer edCo(eds);
  BOOST_CHECK(edCo.getDegree(1)==1);
  BOOST_CHECK(edCo.getDegree(2)==2);
  BOOST_CHECK(edCo.getDegree(3)==1);
  BOOST_CHECK_THROW(edCo.getDegree(4),invalid_argument);
} 

BOOST_AUTO_TEST_CASE(edgeexist_test) {
  Edge ed(1, 2);
  Edge ed2(2, 3);
  vector<Edge> eds{ed, ed2};
  EdgeContainer edCo(eds);
  BOOST_CHECK(edCo.edgeExist(ed));
  BOOST_CHECK(edCo.edgeExist(ed2));
  Edge ed3(3, 4);
  BOOST_CHECK(!edCo.edgeExist(ed3));
}

BOOST_AUTO_TEST_CASE(vertexexist_test) {
  Edge ed(1, 2);
  Edge ed2(2, 3);
  vector<Edge> eds{ed, ed2};
  EdgeContainer edCo(eds);

  BOOST_CHECK(edCo.vertexExist(1));
  BOOST_CHECK(edCo.vertexExist(2));
  BOOST_CHECK(edCo.vertexExist(3));
  BOOST_CHECK(!edCo.vertexExist(4));
}

BOOST_AUTO_TEST_CASE(addedge_test) {
  Edge ed(1, 2);
  Edge ed2(2, 3);
  EdgeContainer edCo;
  edCo.addEdge(ed);
  edCo.addEdge(ed2);
  BOOST_CHECK(edCo.edgeExist(ed));
  BOOST_CHECK(edCo.edgeExist(ed2));
}

BOOST_AUTO_TEST_CASE(getedges_test) {
  Edge ed(1, 2);
  Edge ed2(2, 3);
  EdgeContainer edCo;
  edCo.addEdge(ed);
  edCo.addEdge(ed2);
  auto vec_ed = edCo.getEdges();
  bool ed_found = false;
  bool ed2_found = false;
  for (auto e1 : vec_ed) {
    if (e1 == ed) ed_found = true;
    if (e1 == ed) ed2_found = true;
  }
  BOOST_CHECK(ed_found);
  BOOST_CHECK(ed2_found);
}

BOOST_AUTO_TEST_CASE(getvertices_test) {
  Edge ed(1, 2);
  Edge ed2(2, 3);
  EdgeContainer edCo;
  edCo.addEdge(ed);
  edCo.addEdge(ed2);
  auto vec_vert = edCo.getVertices();
  bool vert_found = false;
  bool vert2_found = false;
  bool vert3_found = false;
  for (auto ver : vec_vert) {
    if (ver == 1) vert_found = true;
    if (ver == 2) vert2_found = true;
    if (ver == 3) vert3_found = true;
  }
  BOOST_CHECK(vert_found);
  BOOST_CHECK(vert2_found);
  BOOST_CHECK(vert3_found);
}

BOOST_AUTO_TEST_CASE(getneighvertices_test) {
  //
  // 1 - 2 - 3
  //
  Edge ed(1, 2);
  Edge ed2(2, 3);
  EdgeContainer edCo;
  edCo.addEdge(ed);
  edCo.addEdge(ed2);
  auto vec_vert = edCo.getNeighVertices(1);
  BOOST_CHECK_EQUAL(vec_vert.at(0), 2);

  vec_vert = edCo.getNeighVertices(2);
  bool vert_found = false;
  bool vert3_found = false;
  for (auto ver : vec_vert) {
    if (ver == 1) vert_found = true;
    if (ver == 3) vert3_found = true;
  }
  BOOST_CHECK(vert_found);
  BOOST_CHECK(vert3_found);
}

BOOST_AUTO_TEST_CASE(getneighedges) {
  Edge ed(1, 2);
  Edge ed2(2, 3);
  EdgeContainer edCo;
  edCo.addEdge(ed);
  edCo.addEdge(ed2);
  auto vec_edgs = edCo.getNeighEdges(1);
  BOOST_CHECK_EQUAL(vec_edgs.at(0), ed);

  vec_edgs = edCo.getNeighEdges(2);
  bool edge_found = false;
  bool edge2_found = false;
  for (auto e1 : vec_edgs) {
    if (e1 == ed) edge_found = true;
    if (e1 == ed2) edge2_found = true;
  }
  BOOST_CHECK(edge_found);
  BOOST_CHECK(edge2_found);
}

BOOST_AUTO_TEST_CASE(getmaxdegree) {
  Edge ed(1,2);
  Edge ed1(2,3);
  Edge ed2(2,4);
  Edge ed3(3,5);
  
  EdgeContainer edCo;
  edCo.addEdge(ed);
  edCo.addEdge(ed1);
  edCo.addEdge(ed2);
  edCo.addEdge(ed3);
  
  int maxD = edCo.getMaxDegree();
  BOOST_CHECK_EQUAL(maxD,3);
}

BOOST_AUTO_TEST_SUITE_END()
