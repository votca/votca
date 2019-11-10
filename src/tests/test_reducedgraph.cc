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

#define BOOST_TEST_MODULE reducedgraph_test
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <exception>
#include <iostream>
#include <votca/tools/graphnode.h>
#include <votca/tools/reducededge.h>
#include <votca/tools/reducedgraph.h>
using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(reducedgraph_test)

BOOST_AUTO_TEST_CASE(constructors_test) { ReducedGraph g; }

/**
 * \brief Test on isolated nodes method
 *
 * The isolated nodes method is meant to grab any nodes that have no edges, as
 * in they exist as islands within the context of the graph.
 */
BOOST_AUTO_TEST_CASE(isolatednodes_test) {

  {

    /// Here gn is a single node as is thus isolated
    vector<ReducedEdge> vec_ed;
    GraphNode gn;
    unordered_map<votca::Index, GraphNode> m_gn;
    m_gn[0] = gn;

    ReducedGraph g(vec_ed, m_gn);
    auto iso_gn = g.getIsolatedNodes();
    BOOST_CHECK_EQUAL(iso_gn.at(0).first, 0);
  }

  {

    /// In this test case gn, gn1 and gn2 are all islands no edges have been
    /// specified to connect them. Calling getIsolatedNodes() thus returns all
    /// three of them.
    vector<ReducedEdge> vec_ed;
    GraphNode gn;
    GraphNode gn1;
    GraphNode gn2;

    unordered_map<votca::Index, GraphNode> m_gn;
    m_gn[0] = gn;
    m_gn[1] = gn1;
    m_gn[2] = gn2;

    ReducedGraph g(vec_ed, m_gn);
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

  {

    /// In this test both node 0 and 1 share an edge and are no longer isolated
    /// however node 2 is isolated, a call getIsolatedNodes() only returns node
    /// 2
    vector<ReducedEdge> vec_ed;
    ReducedEdge ed(0, 1);
    vec_ed.push_back(ed);

    GraphNode gn;
    GraphNode gn1;
    GraphNode gn2;

    unordered_map<votca::Index, GraphNode> m_gn;
    m_gn[0] = gn;
    m_gn[1] = gn1;
    m_gn[2] = gn2;

    ReducedGraph g(vec_ed, m_gn);
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
}

BOOST_AUTO_TEST_CASE(get_edges_test) {

  {
    unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
    unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
    unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
    unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
    unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

    unordered_map<string, double> double_vals;
    unordered_map<string, string> str_vals;

    // 0 - 1 - 2 - 3
    //         |
    //         4

    vector<ReducedEdge> vec_ed;
    ReducedEdge ed(std::vector<votca::Index>{0, 1, 2});
    ReducedEdge ed2(2, 3);
    ReducedEdge ed3(2, 4);

    vec_ed.push_back(ed);
    vec_ed.push_back(ed2);
    vec_ed.push_back(ed3);

    GraphNode gn(int_vals0, double_vals, str_vals);
    GraphNode gn1(int_vals1, double_vals, str_vals);
    GraphNode gn2(int_vals2, double_vals, str_vals);
    GraphNode gn3(int_vals3, double_vals, str_vals);
    GraphNode gn4(int_vals4, double_vals, str_vals);

    unordered_map<votca::Index, GraphNode> m_gn;
    m_gn[0] = gn;
    m_gn[1] = gn1;
    m_gn[2] = gn2;
    m_gn[3] = gn3;
    m_gn[4] = gn4;

    ReducedGraph g(vec_ed, m_gn);
    auto edges = g.getEdges();

    bool ed0_found = false;
    bool ed2_found = false;
    bool ed3_found = false;
    for (auto ed_temp : edges) {
      if (ed_temp == ed) {
        ed0_found = true;
      }
      if (ed_temp == ed2) {
        ed2_found = true;
      }
      if (ed_temp == ed3) {
        ed3_found = true;
      }
    }

    BOOST_CHECK(ed0_found);
    BOOST_CHECK(ed2_found);
    BOOST_CHECK(ed3_found);
  }

  {
    unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
    unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
    unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
    unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
    unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};
    unordered_map<string, votca::Index> int_vals5 = {{"f", 5}};
    unordered_map<string, votca::Index> int_vals6 = {{"g", 6}};
    unordered_map<string, votca::Index> int_vals7 = {{"h", 7}};

    unordered_map<string, double> double_vals;
    unordered_map<string, string> str_vals;

    // Full Graph
    //
    // 0 - 1 - 2 - 3 - 5 - 6
    //     |   |
    //     7 - 4
    //
    //  Reduced Graph
    //
    //  0 - 1 - 2 - 6
    //      | _ |
    //
    vector<ReducedEdge> vec_ed;
    ReducedEdge ed(0, 1);
    ReducedEdge ed2(1, 2);
    ReducedEdge ed3(vector<votca::Index>{2, 3, 4, 5, 6});
    ReducedEdge ed4(vector<votca::Index>{1, 7, 4, 2});

    vec_ed.push_back(ed);
    vec_ed.push_back(ed2);
    vec_ed.push_back(ed3);
    vec_ed.push_back(ed4);

    GraphNode gn(int_vals0, double_vals, str_vals);
    GraphNode gn1(int_vals1, double_vals, str_vals);
    GraphNode gn2(int_vals2, double_vals, str_vals);
    GraphNode gn3(int_vals3, double_vals, str_vals);
    GraphNode gn4(int_vals4, double_vals, str_vals);
    GraphNode gn5(int_vals5, double_vals, str_vals);
    GraphNode gn6(int_vals6, double_vals, str_vals);
    GraphNode gn7(int_vals7, double_vals, str_vals);

    unordered_map<votca::Index, GraphNode> m_gn;
    m_gn[0] = gn;
    m_gn[1] = gn1;
    m_gn[2] = gn2;
    m_gn[3] = gn3;
    m_gn[4] = gn4;
    m_gn[5] = gn5;
    m_gn[6] = gn6;
    m_gn[7] = gn7;

    ReducedGraph g(vec_ed, m_gn);
    auto edges = g.getEdges();

    //  Reduced Graph
    //
    //  0 - 1 - 2 - 6
    //      | _ |
    //
    votca::Index ed0_1count = 0;
    votca::Index ed1_2count = 0;
    votca::Index ed2_6count = 0;

    Edge ed0_1(0, 1);
    Edge ed1_2(1, 2);
    Edge ed2_6(2, 6);

    for (auto ed_temp : edges) {
      if (ed_temp == ed0_1) {
        ++ed0_1count;
      }
      if (ed_temp == ed1_2) {
        ++ed1_2count;
      }
      if (ed_temp == ed2_6) {
        ++ed2_6count;
      }
    }

    BOOST_CHECK_EQUAL(ed0_1count, 1);
    BOOST_CHECK_EQUAL(ed1_2count, 2);
    BOOST_CHECK_EQUAL(ed2_6count, 1);
  }
}

BOOST_AUTO_TEST_CASE(get_vertices_test) {

  unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
  unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
  unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
  unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
  unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

  unordered_map<string, double> double_vals;
  unordered_map<string, string> str_vals;

  // 0 - 1 - 2 - 3
  //         |
  //         4

  vector<ReducedEdge> vec_ed;
  ReducedEdge ed(std::vector<votca::Index>{0, 1, 2});
  ReducedEdge ed2(2, 3);
  ReducedEdge ed3(2, 4);

  vec_ed.push_back(ed);
  vec_ed.push_back(ed2);
  vec_ed.push_back(ed3);

  GraphNode gn(int_vals0, double_vals, str_vals);
  GraphNode gn1(int_vals1, double_vals, str_vals);
  GraphNode gn2(int_vals2, double_vals, str_vals);
  GraphNode gn3(int_vals3, double_vals, str_vals);
  GraphNode gn4(int_vals4, double_vals, str_vals);

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;
  m_gn[3] = gn3;
  m_gn[4] = gn4;

  ReducedGraph g(vec_ed, m_gn);
  vector<votca::Index> vertices = g.getVertices();
  for (auto vertex : vertices) {
    cout << vertex << endl;
  }
  BOOST_CHECK_EQUAL(vertices.size(), 4);
  vector<bool> vertices_found(4, false);
  for (auto vertex : vertices) {
    if (vertex == 0) {
      vertices_found.at(0) = true;
    }
    if (vertex == 2) {
      vertices_found.at(1) = true;
    }
    if (vertex == 3) {
      vertices_found.at(2) = true;
    }
    if (vertex == 4) {
      vertices_found.at(3) = true;
    }
  }

  for (auto found : vertices_found) {
    BOOST_CHECK(found);
  }

  // Lets add a redundant edge, this should not change the output of the
  // getVertices method
  //
  // 0 - 1 - 2 - 3
  //     | _ |
  //         4

  vec_ed.push_back(ed);
  ReducedGraph g2(vec_ed, m_gn);
  vertices = g2.getVertices();
  BOOST_CHECK_EQUAL(vertices.size(), 4);
  vector<bool> vertices_found2(4, false);
  for (auto vertex : vertices) {
    if (vertex == 0) {
      vertices_found2.at(0) = true;
    }
    if (vertex == 2) {
      vertices_found2.at(1) = true;
    }
    if (vertex == 3) {
      vertices_found2.at(2) = true;
    }
    if (vertex == 4) {
      vertices_found2.at(3) = true;
    }
  }

  for (auto found : vertices_found2) {
    BOOST_CHECK(found);
  }
}

BOOST_AUTO_TEST_CASE(compare_test) {
  {

    unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
    unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
    unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
    unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
    unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

    unordered_map<string, double> double_vals;
    unordered_map<string, string> str_vals;

    // 0 - 1 - 2 - 3
    //         |
    //         4

    vector<ReducedEdge> vec_ed;
    ReducedEdge ed(std::vector<votca::Index>{0, 1, 2});
    ReducedEdge ed2(2, 3);
    ReducedEdge ed3(2, 4);

    vec_ed.push_back(ed);
    vec_ed.push_back(ed2);
    vec_ed.push_back(ed3);

    GraphNode gn(int_vals0, double_vals, str_vals);
    GraphNode gn1(int_vals1, double_vals, str_vals);
    GraphNode gn2(int_vals2, double_vals, str_vals);
    GraphNode gn3(int_vals3, double_vals, str_vals);
    GraphNode gn4(int_vals4, double_vals, str_vals);

    unordered_map<votca::Index, GraphNode> m_gn;
    m_gn[0] = gn;
    m_gn[1] = gn1;
    m_gn[2] = gn2;
    m_gn[3] = gn3;
    m_gn[4] = gn4;

    ReducedGraph g(vec_ed, m_gn);

    auto vec_pr = g.getNodes();

    sort(vec_pr.begin(), vec_pr.end(), cmpVertNodePair);
    BOOST_CHECK_EQUAL(vec_pr.at(0).first, 0);
    BOOST_CHECK_EQUAL(vec_pr.at(1).first, 2);
    BOOST_CHECK_EQUAL(vec_pr.at(2).first, 3);
    BOOST_CHECK_EQUAL(vec_pr.at(3).first, 4);
  }

  {
    unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
    unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
    unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
    unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
    unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

    unordered_map<string, double> double_vals;
    unordered_map<string, string> str_vals;

    // 0 - 1 - 2 - 3
    //         |
    //         4
    vector<ReducedEdge> vec_ed;
    ReducedEdge ed(std::vector<votca::Index>{0, 1, 2});
    ReducedEdge ed2(2, 3);
    ReducedEdge ed3(2, 4);

    vec_ed.push_back(ed);
    vec_ed.push_back(ed2);
    vec_ed.push_back(ed3);

    GraphNode gn(int_vals0, double_vals, str_vals);
    GraphNode gn1(int_vals1, double_vals, str_vals);
    GraphNode gn2(int_vals2, double_vals, str_vals);
    GraphNode gn3(int_vals3, double_vals, str_vals);
    GraphNode gn4(int_vals4, double_vals, str_vals);

    /// Only difference is here where we have rearanged the nodes
    unordered_map<votca::Index, GraphNode> m_gn;
    m_gn[4] = gn;
    m_gn[1] = gn1;
    m_gn[3] = gn2;
    m_gn[2] = gn3;
    m_gn[0] = gn4;

    ReducedGraph g(vec_ed, m_gn);

    auto vec_pr = g.getNodes();

    sort(vec_pr.begin(), vec_pr.end(), cmpVertNodePair);
    BOOST_CHECK_EQUAL(vec_pr.at(0).first, 4);
    BOOST_CHECK_EQUAL(vec_pr.at(1).first, 3);
    BOOST_CHECK_EQUAL(vec_pr.at(2).first, 2);
    BOOST_CHECK_EQUAL(vec_pr.at(3).first, 0);
  }
}

BOOST_AUTO_TEST_CASE(neighbornode_test) {
  unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
  unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
  unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
  unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
  unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

  unordered_map<string, double> double_vals;

  unordered_map<string, string> str_vals;

  // 0 - 1 - 2 - 3
  //         |
  //         4
  vector<ReducedEdge> vec_ed;
  ReducedEdge ed(std::vector<votca::Index>{0, 1, 2});
  ReducedEdge ed2(2, 3);
  ReducedEdge ed3(2, 4);

  vec_ed.push_back(ed);
  vec_ed.push_back(ed2);
  vec_ed.push_back(ed3);

  GraphNode gn(int_vals0, double_vals, str_vals);
  GraphNode gn1(int_vals1, double_vals, str_vals);
  GraphNode gn2(int_vals2, double_vals, str_vals);
  GraphNode gn3(int_vals3, double_vals, str_vals);
  GraphNode gn4(int_vals4, double_vals, str_vals);

  unordered_map<votca::Index, GraphNode> m_gn;
  /// Here the graph nodes are assigne to different vertices
  m_gn[0] = gn4;
  m_gn[1] = gn1;
  m_gn[2] = gn3;
  m_gn[3] = gn2;
  m_gn[4] = gn;

  ReducedGraph g(vec_ed, m_gn);

  auto neigh1 = g.getNeighNodes(0);
  BOOST_CHECK_EQUAL(neigh1.size(), 1);
  bool neigh1_found1 = neigh1.at(0).second == gn3;
  BOOST_CHECK(neigh1_found1);

  auto neigh3 = g.getNeighNodes(2);
  BOOST_CHECK_EQUAL(neigh3.size(), 3);
  bool neigh3_found1 = false;
  bool neigh3_found2 = false;
  bool neigh3_found3 = false;
  for (auto neigh_pr : neigh3) {
    if (neigh_pr.second == gn) {
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

BOOST_AUTO_TEST_CASE(expandedge_test) {

  // 0 - 1 - 2 - 3 -4 -5
  //     |       |
  //     6 - 7 - 8
  //
  ReducedEdge ed0(vector<votca::Index>{0, 1});
  ReducedEdge ed1(vector<votca::Index>{1, 2, 3});
  ReducedEdge ed2(vector<votca::Index>{3, 4, 5});
  ReducedEdge ed3(vector<votca::Index>{1, 6, 7, 8, 3});

  vector<ReducedEdge> vec_ed{ed0, ed1, ed2, ed3};

  ReducedGraph g(vec_ed);

  vector<vector<Edge>> edges = g.expandEdge(Edge(0, 1));
  BOOST_CHECK_EQUAL(edges.size(), 1);
  BOOST_CHECK_EQUAL(edges.at(0).size(), 1);

  edges = g.expandEdge(Edge(1, 3));
  BOOST_CHECK_EQUAL(edges.size(), 2);  // Two separate chains

  if (edges.at(0).size() == 2) {
    bool found_ed1_2 = false;
    bool found_ed2_3 = false;

    Edge ed1_2(1, 2);
    Edge ed2_3(2, 3);
    for (auto e1 : edges.at(0)) {
      if (e1 == ed1_2) {
        found_ed1_2 = true;
      }
      if (e1 == ed2_3) {
        found_ed2_3 = true;
      }
    }
    BOOST_CHECK(found_ed1_2);
    BOOST_CHECK(found_ed2_3);

    Edge ed1_6(1, 6);
    Edge ed6_7(6, 7);
    Edge ed7_8(7, 8);
    Edge ed8_3(3, 8);

    bool found_ed1_6 = false;
    bool found_ed6_7 = false;
    bool found_ed7_8 = false;
    bool found_ed8_3 = false;

    for (auto e1 : edges.at(1)) {
      if (e1 == ed1_6) {
        found_ed1_6 = true;
      }
      if (e1 == ed6_7) {
        found_ed6_7 = true;
      }
      if (e1 == ed7_8) {
        found_ed7_8 = true;
      }
      if (e1 == ed8_3) {
        found_ed8_3 = true;
      }
    }
    BOOST_CHECK(found_ed1_6);
    BOOST_CHECK(found_ed6_7);
    BOOST_CHECK(found_ed7_8);
    BOOST_CHECK(found_ed8_3);

  } else if (edges.at(0).size() == 4) {
    bool found_ed1_2 = false;
    bool found_ed2_3 = false;

    Edge ed1_2(1, 2);
    Edge ed2_3(2, 3);
    for (auto e1 : edges.at(1)) {
      if (e1 == ed1_2) {
        found_ed1_2 = true;
      }
      if (e1 == ed2_3) {
        found_ed2_3 = true;
      }
    }
    BOOST_CHECK(found_ed1_2);
    BOOST_CHECK(found_ed2_3);

    Edge ed1_6(1, 6);
    Edge ed6_7(6, 7);
    Edge ed7_8(7, 8);
    Edge ed8_3(3, 8);

    bool found_ed1_6 = false;
    bool found_ed6_7 = false;
    bool found_ed7_8 = false;
    bool found_ed8_3 = false;

    for (auto e1 : edges.at(0)) {
      if (e1 == ed1_6) {
        found_ed1_6 = true;
      }
      if (e1 == ed6_7) {
        found_ed6_7 = true;
      }
      if (e1 == ed7_8) {
        found_ed7_8 = true;
      }
      if (e1 == ed8_3) {
        found_ed8_3 = true;
      }
    }
    BOOST_CHECK(found_ed1_6);
    BOOST_CHECK(found_ed6_7);
    BOOST_CHECK(found_ed7_8);
    BOOST_CHECK(found_ed8_3);

  } else {
    // ONe of the two options above should have been triggered
    BOOST_CHECK(false);
  }
}

BOOST_AUTO_TEST_CASE(getdegree_test) {

  //             9
  //             |
  // 0 - 1 - 2 - 3 -4 -5
  //     |       |
  //     6 - 7 - 8
  //
  ReducedEdge ed0(vector<votca::Index>{0, 1});
  ReducedEdge ed1(vector<votca::Index>{1, 2, 3});
  ReducedEdge ed2(vector<votca::Index>{3, 4, 5});
  ReducedEdge ed3(vector<votca::Index>{1, 6, 7, 8, 3});
  ReducedEdge ed4(3, 9);

  vector<ReducedEdge> vec_ed{ed0, ed1, ed2, ed3, ed4};

  ReducedGraph g(vec_ed);

  BOOST_CHECK_EQUAL(g.getMaxDegree(), 4);

  BOOST_CHECK_EQUAL(g.getDegree(0), 1);
  BOOST_CHECK_EQUAL(g.getDegree(1), 3);
  BOOST_CHECK_EQUAL(g.getDegree(3), 4);
  BOOST_CHECK_EQUAL(g.getDegree(5), 1);
  BOOST_CHECK_EQUAL(g.getDegree(9), 1);

  votca::Index degree = 0;
  auto vertices_w_degree = g.getVerticesDegree(degree);
  BOOST_CHECK_EQUAL(vertices_w_degree.size(), 0);

  degree = 1;
  vertices_w_degree = g.getVerticesDegree(degree);
  BOOST_CHECK_EQUAL(vertices_w_degree.size(), 3);

  degree = 2;
  vertices_w_degree = g.getVerticesDegree(degree);
  BOOST_CHECK_EQUAL(vertices_w_degree.size(), 0);

  degree = 3;
  vertices_w_degree = g.getVerticesDegree(degree);
  BOOST_CHECK_EQUAL(vertices_w_degree.size(), 1);

  degree = 4;
  vertices_w_degree = g.getVerticesDegree(degree);
  BOOST_CHECK_EQUAL(vertices_w_degree.size(), 1);
}

/**
 * \brief Equivalence test
 *
 * Here we demonstrate how the equivalence test works it is purely dependendent
 * on whether the contents of the graphnodes in the graph contain the same
 * information.
 */
BOOST_AUTO_TEST_CASE(id_test) {
  {
    unordered_map<string, votca::Index> int_vals0 = {{"a", 0}};
    unordered_map<string, votca::Index> int_vals1 = {{"b", 1}};
    unordered_map<string, votca::Index> int_vals2 = {{"c", 2}};
    unordered_map<string, votca::Index> int_vals3 = {{"d", 3}};
    unordered_map<string, votca::Index> int_vals4 = {{"e", 4}};

    unordered_map<string, double> double_vals;

    unordered_map<string, string> str_vals;

    vector<ReducedEdge> vec_ed;
    ReducedEdge ed(std::vector<votca::Index>{0, 1, 2});
    ReducedEdge ed2(2, 3);
    ReducedEdge ed3(2, 4);

    vec_ed.push_back(ed);
    vec_ed.push_back(ed2);
    vec_ed.push_back(ed3);

    GraphNode gn(int_vals0, double_vals, str_vals);
    GraphNode gn1(int_vals1, double_vals, str_vals);
    GraphNode gn2(int_vals2, double_vals, str_vals);
    GraphNode gn3(int_vals3, double_vals, str_vals);
    GraphNode gn4(int_vals4, double_vals, str_vals);

    unordered_map<votca::Index, GraphNode> m_gn;
    /// Here the graph nodes are assigne to different vertices
    m_gn[4] = gn;
    m_gn[1] = gn1;
    m_gn[3] = gn2;
    m_gn[2] = gn3;
    m_gn[0] = gn4;

    ReducedGraph g(vec_ed, m_gn);

    /// Here is what the string id of the graph should look like
    string str = "a0c2d3e4";
    string s_id = g.getId();
    BOOST_CHECK_EQUAL(s_id, str);

    ReducedGraph g2(vec_ed, m_gn);
    BOOST_CHECK(g == g2);

    /// Here we switch up which vertices contain which graphnodes and show that
    /// the graph id is the same. This is because the vertex ids are not used to
    /// create the id and neither the edges. Only the contens in the graphnodes
    m_gn[1] = gn3;
    m_gn[2] = gn1;
    ReducedGraph g3(vec_ed, m_gn);
    BOOST_CHECK(g != g3);

    string str2 = "a0b1c2e4";
    string s_id2 = g3.getId();
    BOOST_CHECK_EQUAL(s_id2, str2);

    GraphNode gn5(int_vals3, double_vals, str_vals);
    m_gn[5] = gn5;
    ReducedGraph g4(vec_ed, m_gn);
    BOOST_CHECK(g != g4);
  }
}

BOOST_AUTO_TEST_SUITE_END()
