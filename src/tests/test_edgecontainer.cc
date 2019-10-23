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

#define BOOST_TEST_MODULE edgecontainer_test
#include <boost/test/unit_test.hpp>
#include <exception>
#include <iostream>
#include <vector>
#include <votca/tools/edge.h>
#include <votca/tools/edgecontainer.h>

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

BOOST_AUTO_TEST_CASE(getdegree_test) {
  {
    Edge ed(1, 2);
    Edge ed2(2, 3);
    vector<Edge> eds{ed, ed2};
    EdgeContainer edge_container(eds);
    BOOST_CHECK(edge_container.getDegree(1) == 1);
    BOOST_CHECK(edge_container.getDegree(2) == 2);
    BOOST_CHECK(edge_container.getDegree(3) == 1);
    BOOST_CHECK_THROW(edge_container.getDegree(4), invalid_argument);
  }

  {
    //    _
    //  /  \      3
    //  \  /      |
    //    1 - 4 - 8 - 5
    //
    Edge ed1(1, 4);
    Edge ed2(1, 1);
    Edge ed3(4, 8);
    Edge ed4(8, 3);
    Edge ed5(8, 5);
    vector<Edge> edges{ed1, ed2, ed3, ed4, ed5};
    EdgeContainer edge_container(edges);

    BOOST_CHECK_EQUAL(edge_container.getDegree(1), 3);
    BOOST_CHECK_EQUAL(edge_container.getDegree(4), 2);
    BOOST_CHECK_EQUAL(edge_container.getDegree(8), 3);
    BOOST_CHECK_EQUAL(edge_container.getDegree(3), 1);
    BOOST_CHECK_EQUAL(edge_container.getDegree(5), 1);
  }
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

  // Should be able to add the same edge more than once
  edCo.addEdge(ed);
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
    if (e1 == ed) {
      ed_found = true;
    }
    if (e1 == ed2) {
      ed2_found = true;
    }
  }
  BOOST_CHECK(ed_found);
  BOOST_CHECK(ed2_found);

  // Should be able to add an edge more than once
  edCo.addEdge(ed);
  vec_ed = edCo.getEdges();
  int ed_count = 0;
  int ed2_count = 0;
  for (auto e1 : vec_ed) {
    if (e1 == ed) {
      ++ed_count;
    }
    if (e1 == ed2) {
      ++ed2_count;
    }
  }
  BOOST_CHECK_EQUAL(ed_count, 2);
  BOOST_CHECK_EQUAL(ed2_count, 1);
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
    if (ver == 1) {
      vert_found = true;
    }
    if (ver == 2) {
      vert2_found = true;
    }
    if (ver == 3) {
      vert3_found = true;
    }
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
    if (ver == 1) {
      vert_found = true;
    }
    if (ver == 3) {
      vert3_found = true;
    }
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
    if (e1 == ed) {
      edge_found = true;
    }
    if (e1 == ed2) {
      edge2_found = true;
    }
  }
  BOOST_CHECK(edge_found);
  BOOST_CHECK(edge2_found);

  // Should be able to add the same edge more than once
  Edge ed3(3, 4);
  edCo.addEdge(ed);
  edCo.addEdge(ed3);

  vec_edgs = edCo.getNeighEdges(1);
  BOOST_CHECK_EQUAL(vec_edgs.size(), 2);
  BOOST_CHECK_EQUAL(vec_edgs.at(0), ed);
  BOOST_CHECK_EQUAL(vec_edgs.at(1), ed);

  vec_edgs = edCo.getNeighEdges(2);
  BOOST_CHECK_EQUAL(vec_edgs.size(), 3);

  int edge_count = 0;
  int edge_count2 = 0;
  for (auto e1 : vec_edgs) {
    if (e1 == ed) {
      ++edge_count;
    }
    if (e1 == ed2) {
      ++edge_count2;
    }
  }
  BOOST_CHECK_EQUAL(edge_count, 2);
  BOOST_CHECK_EQUAL(edge_count2, 1);
}

BOOST_AUTO_TEST_CASE(getmaxdegree) {
  Edge ed(1, 2);
  Edge ed1(2, 3);
  Edge ed2(2, 4);
  Edge ed3(3, 5);

  EdgeContainer edCo;
  edCo.addEdge(ed);
  edCo.addEdge(ed1);
  edCo.addEdge(ed2);
  edCo.addEdge(ed3);

  long int maxD = edCo.getMaxDegree();
  BOOST_CHECK_EQUAL(maxD, 3);

  edCo.addEdge(ed);
  maxD = edCo.getMaxDegree();
  BOOST_CHECK_EQUAL(maxD, 4);
}

BOOST_AUTO_TEST_SUITE_END()
