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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE graph_test

// Standard includes
#include <cmath>
#include <exception>
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/edge.h"
#include "votca/tools/graph.h"
#include "votca/tools/graphnode.h"

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(graph_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Graph g; }

/**
 * \brief Test on isolated nodes method
 *
 * The isolated nodes method is meant to grab any nodes that have no edges, as
 * in they exist as islands within the context of the graph.
 */
BOOST_AUTO_TEST_CASE(isolatednodes_test1) {
  /// Here gn is a single node as is thus isolated
  //
  //  gn
  //
  vector<Edge> vec_ed;
  GraphNode gn;
  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;

  Graph g(vec_ed, m_gn);
  auto iso_gn = g.getIsolatedNodes();
  BOOST_CHECK_EQUAL(iso_gn.at(0).first, 0);
}

BOOST_AUTO_TEST_CASE(isolatednodes_test2) {
  /// In this test case gn, gn1 and gn2 are all islands no edges have been
  /// specified to connect them. Calling getIsolatedNodes() thus returns all
  /// three of them.
  //
  //   gn   gn1   gn2
  //
  vector<Edge> vec_ed;
  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;

  Graph g(vec_ed, m_gn);
  auto iso_gn = g.getIsolatedNodes();
  bool node0 = false;
  bool node1 = false;
  bool node2 = false;

  for (auto n_pr : iso_gn) {
    if (n_pr.first == 0) {
      node0 = true;
    }
    if (n_pr.first == 1) {
      node1 = true;
    }
    if (n_pr.first == 2) {
      node2 = true;
    }
  }

  BOOST_CHECK(node0);
  BOOST_CHECK(node1);
  BOOST_CHECK(node2);
}

BOOST_AUTO_TEST_CASE(isolatednodes_test3) {

  /// In this test both node 0 and 1 share an edge and are no longer isolated
  /// however node 2 is isolated, a call getIsolatedNodes() only returns node
  /// 2
  //
  //   gn - - gn1    gn2
  //
  vector<Edge> vec_ed;
  Edge ed(0, 1);
  vec_ed.push_back(ed);

  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;

  Graph g(vec_ed, m_gn);
  auto iso_gn = g.getIsolatedNodes();
  bool node0 = false;
  bool node1 = false;
  bool node2 = false;

  for (auto n_pr : iso_gn) {
    if (n_pr.first == 0) {
      node0 = true;
    }
    if (n_pr.first == 1) {
      node1 = true;
    }
    if (n_pr.first == 2) {
      node2 = true;
    }
  }

  BOOST_CHECK(!node0);
  BOOST_CHECK(!node1);
  BOOST_CHECK(node2);
}

BOOST_AUTO_TEST_CASE(get_junctions_test) {
  /// In this test the junctions of the graph are returned, the junctions
  /// consist of vertices of three or more connections
  /// In this test both node 0 and 1 share an edge there are no junctions
  ///
  /// gn - - gn1   gn2
  ///
  Edge ed(0, 1);
  vector<Edge> vec_ed{ed};

  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;

  Graph g(vec_ed, m_gn);
  auto junctions = g.getJunctions();

  BOOST_CHECK_EQUAL(junctions.size(), 0);
}

BOOST_AUTO_TEST_CASE(get_junctions_test2) {
  /// In this test the junctions of the graph are returned, the junctions
  /// consist of vertices of three or more connections

  /// In this test both node 0, 1 and 3 are all connected to node 2 
  ///
  ///  gn - - gn2 - - gn3
  ///          |
  ///         gn1
  ///
  Edge ed(0, 2);
  Edge ed2(1, 2);
  Edge ed3(3, 2);
  vector<Edge> vec_ed{ed, ed2, ed3};

  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;
  GraphNode gn3;

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;
  m_gn[3] = gn3;

  Graph g(vec_ed, m_gn);
  auto junctions = g.getJunctions();

  BOOST_CHECK_EQUAL(junctions.size(), 1);
  BOOST_CHECK_EQUAL(junctions.at(0), 2);
}

BOOST_AUTO_TEST_CASE(junctions_test3) {

  unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
  unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
  unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
  unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
  unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

  // 0 - 1 - 2 - 3
  //         |
  //         4

  vector<Edge> vec_ed;
  Edge ed(0, 1);
  Edge ed1(1, 2);
  Edge ed2(2, 3);
  Edge ed3(2, 4);

  vec_ed.push_back(ed);
  vec_ed.push_back(ed1);
  vec_ed.push_back(ed2);
  vec_ed.push_back(ed3);

  GraphNode gn(int_vals0);
  GraphNode gn1(int_vals1);
  GraphNode gn2(int_vals2);
  GraphNode gn3(int_vals3);
  GraphNode gn4(int_vals4);

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;
  m_gn[3] = gn3;
  m_gn[4] = gn4;

  Graph g(vec_ed, m_gn);
  // A junction should consist of a vertex with degree of 3 or more
  auto junctions = g.getJunctions();
  BOOST_CHECK(junctions.size() == 1);
  BOOST_CHECK_EQUAL(junctions.at(0), 2);
}

BOOST_AUTO_TEST_CASE(get_edges_test) {

  unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
  unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
  unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
  unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
  unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

  // 0 - 1 - 2 - 3
  //         |
  //         4

  vector<Edge> vec_ed;
  Edge ed(0, 1);
  Edge ed1(1, 2);
  Edge ed2(2, 3);
  Edge ed3(2, 4);

  vec_ed.push_back(ed);
  vec_ed.push_back(ed1);
  vec_ed.push_back(ed2);
  vec_ed.push_back(ed3);

  GraphNode gn(int_vals0);
  GraphNode gn1(int_vals1);
  GraphNode gn2(int_vals2);
  GraphNode gn3(int_vals3);
  GraphNode gn4(int_vals4);

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;
  m_gn[3] = gn3;
  m_gn[4] = gn4;

  Graph g(vec_ed, m_gn);
  auto edges = g.getEdges();

  bool ed0_found = false;
  bool ed1_found = false;
  bool ed2_found = false;
  bool ed3_found = false;
  for (auto ed_temp : edges) {
    if (ed_temp == ed) {
      ed0_found = true;
    }
    if (ed_temp == ed1) {
      ed1_found = true;
    }
    if (ed_temp == ed2) {
      ed2_found = true;
    }
    if (ed_temp == ed3) {
      ed3_found = true;
    }
  }

  BOOST_CHECK(ed0_found);
  BOOST_CHECK(ed1_found);
  BOOST_CHECK(ed2_found);
  BOOST_CHECK(ed3_found);
}

BOOST_AUTO_TEST_CASE(get_vertices_test) {

  unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
  unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
  unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
  unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
  unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

  // 0 - 1 - 2 - 3
  //         |
  //         4

  vector<Edge> vec_ed;
  Edge ed(0, 1);
  Edge ed1(1, 2);
  Edge ed2(2, 3);
  Edge ed3(2, 4);

  vec_ed.push_back(ed);
  vec_ed.push_back(ed1);
  vec_ed.push_back(ed2);
  vec_ed.push_back(ed3);

  GraphNode gn(int_vals0);
  GraphNode gn1(int_vals1);
  GraphNode gn2(int_vals2);
  GraphNode gn3(int_vals3);
  GraphNode gn4(int_vals4);

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;
  m_gn[3] = gn3;
  m_gn[4] = gn4;

  Graph g(vec_ed, m_gn);
  auto vertices = g.getVertices();

  vector<bool> vertices_found(5, false);
  for (auto vertex : vertices) {
    if (vertex == 0) {
      vertices_found.at(0) = true;
    }
    if (vertex == 1) {
      vertices_found.at(1) = true;
    }
    if (vertex == 2) {
      vertices_found.at(2) = true;
    }
    if (vertex == 3) {
      vertices_found.at(3) = true;
    }
    if (vertex == 4) {
      vertices_found.at(4) = true;
    }
  }

  for (auto found : vertices_found) {
    BOOST_CHECK(found);
  }
}

BOOST_AUTO_TEST_CASE(compare_test1) {

  unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
  unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
  unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
  unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
  unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

  // 0 - 1 - 2 - 3
  //         |
  //         4

  vector<Edge> vec_ed;
  Edge ed(0, 1);
  Edge ed1(1, 2);
  Edge ed2(2, 3);
  Edge ed3(2, 4);

  vec_ed.push_back(ed);
  vec_ed.push_back(ed1);
  vec_ed.push_back(ed2);
  vec_ed.push_back(ed3);

  GraphNode gn(int_vals0);
  GraphNode gn1(int_vals1);
  GraphNode gn2(int_vals2);
  GraphNode gn3(int_vals3);
  GraphNode gn4(int_vals4);

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;
  m_gn[3] = gn3;
  m_gn[4] = gn4;

  Graph g(vec_ed, m_gn);

  auto vec_pr = g.getNodes();

  sort(vec_pr.begin(), vec_pr.end(), cmpVertNodePair);
  BOOST_CHECK_EQUAL(vec_pr.at(0).first, 0);
  BOOST_CHECK_EQUAL(vec_pr.at(1).first, 1);
  BOOST_CHECK_EQUAL(vec_pr.at(2).first, 2);
  BOOST_CHECK_EQUAL(vec_pr.at(3).first, 3);
  BOOST_CHECK_EQUAL(vec_pr.at(4).first, 4);
}

BOOST_AUTO_TEST_CASE(compare_test2) {
  unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
  unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
  unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
  unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
  unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

  // 0 - 1 - 2 - 3
  //         |
  //         4
  vector<Edge> vec_ed;
  Edge ed(0, 1);
  Edge ed1(1, 2);
  Edge ed2(2, 3);
  Edge ed3(2, 4);

  vec_ed.push_back(ed);
  vec_ed.push_back(ed1);
  vec_ed.push_back(ed2);
  vec_ed.push_back(ed3);

  GraphNode gn(int_vals0);
  GraphNode gn1(int_vals1);
  GraphNode gn2(int_vals2);
  GraphNode gn3(int_vals3);
  GraphNode gn4(int_vals4);

  /// Only difference is here where we have rearanged the nodes
  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[4] = gn;
  m_gn[1] = gn1;
  m_gn[3] = gn2;
  m_gn[2] = gn3;
  m_gn[0] = gn4;

  Graph g(vec_ed, m_gn);

  auto vec_pr = g.getNodes();

  sort(vec_pr.begin(), vec_pr.end(), cmpVertNodePair);
  BOOST_CHECK_EQUAL(vec_pr.at(0).first, 4);
  BOOST_CHECK_EQUAL(vec_pr.at(1).first, 1);
  BOOST_CHECK_EQUAL(vec_pr.at(2).first, 3);
  BOOST_CHECK_EQUAL(vec_pr.at(3).first, 2);
  BOOST_CHECK_EQUAL(vec_pr.at(4).first, 0);
}

BOOST_AUTO_TEST_CASE(neighbornode_test) {
  unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
  unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
  unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
  unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
  unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

  // 0 - 1 - 2 - 3
  //         |
  //         4
  vector<Edge> vec_ed;
  Edge ed(0, 1);
  Edge ed1(1, 2);
  Edge ed2(2, 3);
  Edge ed3(2, 4);

  vec_ed.push_back(ed);
  vec_ed.push_back(ed1);
  vec_ed.push_back(ed2);
  vec_ed.push_back(ed3);

  GraphNode gn(int_vals0);
  GraphNode gn1(int_vals1);
  GraphNode gn2(int_vals2);
  GraphNode gn3(int_vals3);
  GraphNode gn4(int_vals4);

  unordered_map<votca::Index, GraphNode> m_gn;
  /// Here the graph nodes are assigne to different vertices
  m_gn[0] = gn4;
  m_gn[1] = gn1;
  m_gn[2] = gn3;
  m_gn[3] = gn2;
  m_gn[4] = gn;

  Graph g(vec_ed, m_gn);

  auto neigh1 = g.getNeighNodes(0);
  BOOST_CHECK_EQUAL(neigh1.size(), 1);
  bool neigh1_found1 = neigh1.at(0).second == gn1;
  BOOST_CHECK(neigh1_found1);

  auto neigh2 = g.getNeighNodes(1);
  BOOST_CHECK_EQUAL(neigh2.size(), 2);
  bool neigh2_found1 = false;
  bool neigh2_found2 = false;
  for (auto neigh_pr : neigh2) {
    if (neigh_pr.second == gn4) {
      neigh2_found1 = true;
    }
    if (neigh_pr.second == gn3) {
      neigh2_found2 = true;
    }
  }
  BOOST_CHECK(neigh2_found1);
  BOOST_CHECK(neigh2_found2);

  auto neigh3 = g.getNeighNodes(2);
  BOOST_CHECK_EQUAL(neigh3.size(), 3);
  bool neigh3_found1 = false;
  bool neigh3_found2 = false;
  bool neigh3_found3 = false;
  for (auto neigh_pr : neigh3) {
    if (neigh_pr.second == gn1) {
      neigh3_found1 = true;
    }
    if (neigh_pr.second == gn2) {
      neigh3_found2 = true;
    }
    if (neigh_pr.second == gn) {
      neigh3_found3 = true;
    }
  }
  BOOST_CHECK(neigh3_found1);
  BOOST_CHECK(neigh3_found2);
  BOOST_CHECK(neigh3_found3);

  auto neigh4 = g.getNeighNodes(3);
  BOOST_CHECK_EQUAL(neigh4.size(), 1);
  bool neigh4_found1 = neigh4.at(0).second == gn3;
  BOOST_CHECK(neigh4_found1);

  auto neigh5 = g.getNeighNodes(4);
  BOOST_CHECK_EQUAL(neigh5.size(), 1);
  bool neigh5_found1 = neigh5.at(0).second == gn3;
  BOOST_CHECK(neigh5_found1);
}

/**
 * \brief Equivalence test
 *
 * Here we demonstrate how the equivalence test works it is purely dependendent
 * on whether the contents of the graphnodes in the graph contain the same
 * information.
 */
BOOST_AUTO_TEST_CASE(id_test) {
  unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
  unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
  unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
  unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
  unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

  vector<Edge> vec_ed;
  Edge ed(0, 1);
  Edge ed1(1, 2);
  Edge ed2(2, 3);
  Edge ed3(2, 4);

  vec_ed.push_back(ed);
  vec_ed.push_back(ed1);
  vec_ed.push_back(ed2);
  vec_ed.push_back(ed3);

  GraphNode gn(int_vals0);
  GraphNode gn1(int_vals1);
  GraphNode gn2(int_vals2);
  GraphNode gn3(int_vals3);
  GraphNode gn4(int_vals4);

  unordered_map<votca::Index, GraphNode> m_gn;
  /// Here the graph nodes are assigne to different vertices
  m_gn[4] = gn;
  m_gn[1] = gn1;
  m_gn[3] = gn2;
  m_gn[2] = gn3;
  m_gn[0] = gn4;

  Graph g(vec_ed, m_gn);

  /// Here is what the string id of the graph should look like
  string str = "a=0;b=1;c=2;d=3;e=4;";
  string s_id = g.getContentLabel().get();
  BOOST_CHECK_EQUAL(s_id, str);

  Graph g2(vec_ed, m_gn);
  BOOST_CHECK(g == g2);

  /// Here we switch up which vertices contain which graphnodes and show that
  /// the graph id is the same. This is because the vertex ids are not used to
  /// create the id and neither the edges. Only the contens in the graphnodes
  m_gn[1] = gn3;
  m_gn[2] = gn1;
  Graph g3(vec_ed, m_gn);
  BOOST_CHECK(g == g3);

  GraphNode gn5(int_vals3);
  m_gn[5] = gn5;
  Graph g4(vec_ed, m_gn);
  BOOST_CHECK(g != g4);
}

BOOST_AUTO_TEST_SUITE_END()
