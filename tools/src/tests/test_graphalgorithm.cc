/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE graphalgorithm_test

// Standard includes
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/edge.h"
#include "votca/tools/graph.h"
#include "votca/tools/graph_bf_visitor.h"
#include "votca/tools/graphalgorithm.h"
#include "votca/tools/graphdistvisitor.h"
#include "votca/tools/graphnode.h"
#include "votca/tools/reducedgraph.h"

using namespace std;
using namespace votca::tools;
using namespace boost;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(graphalgorithm_test)

BOOST_AUTO_TEST_CASE(single_network_algorithm_test) {
  /**
   * @brief The singleNetwork funciton is used to determine if a graph is 
   * composed of a single connected network. To satisfy this criteria all
   * of the vertices in the graph must be connected through a finite number
   * of connections.
   *
   * In this test we add two nodes and an edge describing
   * their connection thus the singleNetwork function will
   * return true.
   *
   *   gn1 - - gn2
   *
   */

  // Create edge
  Edge ed(0, 1);
  vector<Edge> edges;
  edges.push_back(ed);

  // Create Graph nodes
  GraphNode gn1;
  GraphNode gn2;

  unordered_map<votca::Index, GraphNode> nodes;
  nodes[0] = gn1;
  nodes[1] = gn2;

  Graph g(edges, nodes);

  Graph_BF_Visitor gb_v;

  BOOST_CHECK(gb_v.queEmpty());
  BOOST_CHECK_THROW(gb_v.exec(&g, ed), runtime_error);

  bool single_n = singleNetwork(g, gb_v);
  cerr << "is single network " << single_n << endl;
  BOOST_CHECK(single_n);
  BOOST_CHECK(gb_v.queEmpty());

  Graph_BF_Visitor gb_v2;
  gb_v2.setStartingVertex(2);
  BOOST_CHECK_THROW(singleNetwork(g, gb_v2), invalid_argument);
}

BOOST_AUTO_TEST_CASE(single_network_algorithm_test2) {

  /** In this test we add 3 nodes but only one edge
   * this means that one of the nodes will not be
   * attached to the other two. Thus the singleNetwork
   * function should return false.
   *
   * gn1 - - gn2     gn3
   */

  // Create edge
  Edge ed(0, 1);
  vector<Edge> edges;
  edges.push_back(ed);

  // Create Graph nodes
  GraphNode gn1;
  GraphNode gn2;
  GraphNode gn3;

  unordered_map<votca::Index, GraphNode> nodes;
  nodes[0] = gn1;
  nodes[1] = gn2;
  nodes[2] = gn3;

  Graph g(edges, nodes);

  Graph_BF_Visitor gb_v;

  BOOST_CHECK(gb_v.queEmpty());
  BOOST_CHECK_THROW(gb_v.exec(&g, ed), runtime_error);

  bool single_n = singleNetwork(g, gb_v);
  BOOST_CHECK(!single_n);
  BOOST_CHECK(gb_v.queEmpty());
}

BOOST_AUTO_TEST_CASE(decouple_isolated_subgraphs_test) {

  /*
   * The decoupleIsolatedSubGraphs functions is designed to split up a graph
   * based on whether the networks in the graph are connected. 
   *
   * E.g. given a single graph consisting of the following vertices and 
   * edges. 
   *
   *  1 - 2 - 3
   *      |   |            8 - 9 - 10      11
   *      4 - 5 - 6 -7
   *
   * Calling `decoupleIsolatedSubGraphs` on this graph will result in the
   * generation of three separate subgraphs
   */
  Edge ed1(1, 2);
  Edge ed2(2, 3);
  Edge ed3(2, 4);
  Edge ed4(4, 5);
  Edge ed5(3, 5);
  Edge ed6(5, 6);
  Edge ed7(6, 7);

  Edge ed8(8, 9);
  Edge ed9(9, 10);

  vector<Edge> edges{ed1, ed2, ed3, ed4, ed5, ed6, ed7, ed8, ed9};

  unordered_map<votca::Index, GraphNode> nodes;
  for (votca::Index index = 1; index < 12; ++index) {
    GraphNode gn;
    nodes[index] = gn;
  }

  Graph graph(edges, nodes);

  vector<Graph> sub_graphs = votca::tools::decoupleIsolatedSubGraphs(graph);

  BOOST_CHECK_EQUAL(sub_graphs.size(), 3);

  unordered_map<votca::Index, bool> vertices_sub_graph1;
  vertices_sub_graph1[1] = false;
  vertices_sub_graph1[2] = false;
  vertices_sub_graph1[3] = false;
  vertices_sub_graph1[4] = false;
  vertices_sub_graph1[5] = false;
  vertices_sub_graph1[6] = false;
  vertices_sub_graph1[7] = false;

  unordered_map<votca::Index, bool> vertices_sub_graph2;
  vertices_sub_graph2[8] = false;
  vertices_sub_graph2[9] = false;
  vertices_sub_graph2[10] = false;
  unordered_map<votca::Index, bool> vertices_sub_graph3;
  vertices_sub_graph3[11] = false;

  for (const Graph& sub_graph : sub_graphs) {
    const vector<votca::Index> vertices = sub_graph.getVertices();
    if (vertices_sub_graph1.count(vertices.at(0))) {
      /* Checking sub graph 1
       *
       *  1 - 2 - 3
       *      |   |
       *      4 - 5 - 6 -7
       */
      for (const votca::Index& vertex : vertices) {
        if (vertices_sub_graph1.count(vertex)) {
          vertices_sub_graph1.at(vertex) = true;
        }
      }
    } else if (vertices_sub_graph2.count(vertices.at(0))) {
      /* Checking sub graph 2
       *
       *  8 - 9 - 10
       */
      for (const votca::Index& vertex : vertices) {
        if (vertices_sub_graph2.count(vertex)) {
          vertices_sub_graph2.at(vertex) = true;
        }
      }
    } else if (vertices_sub_graph3.count(vertices.at(0))) {
      /* Checking sub graph 3
       * 
       *   11
       */
      for (const votca::Index& vertex : vertices) {
        if (vertices_sub_graph3.count(vertex)) {
          vertices_sub_graph3.at(vertex) = true;
        }
      }
    }
  }

  BOOST_CHECK_EQUAL(vertices_sub_graph1.size(), 7);
  for (const pair<const votca::Index, bool>& found : vertices_sub_graph1) {
    BOOST_CHECK(found.second);
  }

  BOOST_CHECK_EQUAL(vertices_sub_graph2.size(), 3);
  for (const pair<const votca::Index, bool>& found : vertices_sub_graph2) {
    BOOST_CHECK(found.second);
  }

  BOOST_CHECK_EQUAL(vertices_sub_graph3.size(), 1);
  for (const pair<const votca::Index, bool>& found : vertices_sub_graph3) {
    BOOST_CHECK(found.second);
  }
}

BOOST_AUTO_TEST_CASE(reduceGraph_algorithm_test) {

  // Create edge
  //     2
  //     |
  // 0 - 1 - 3
  //     |
  //     4

  Edge ed1(0, 1);
  Edge ed2(1, 2);
  Edge ed3(1, 3);
  Edge ed4(1, 4);

  //
  // 5 - 6
  //
  Edge ed5(5, 6);

  vector<Edge> edges;
  edges.push_back(ed1);
  edges.push_back(ed2);
  edges.push_back(ed3);
  edges.push_back(ed4);
  edges.push_back(ed5);

  // Create Graph nodes
  GraphNode gn1;
  GraphNode gn2;
  GraphNode gn3;
  GraphNode gn4;
  GraphNode gn5;
  GraphNode gn6;
  GraphNode gn7;

  unordered_map<votca::Index, GraphNode> nodes;
  nodes[0] = gn1;
  nodes[1] = gn2;
  nodes[2] = gn3;
  nodes[3] = gn4;
  nodes[4] = gn5;
  nodes[5] = gn6;
  nodes[6] = gn7;

  Graph g(edges, nodes);
  ReducedGraph reduced_g = votca::tools::reduceGraph(g);

  vector<Edge> edges2 = reduced_g.getEdges();
  BOOST_CHECK_EQUAL(edges2.size(), 5);
  vector<bool> found_edges(5, false);
  for (Edge& edge : edges2) {
    if (edge == ed1) {
      found_edges.at(0) = true;
    }
    if (edge == ed2) {
      found_edges.at(1) = true;
    }
    if (edge == ed3) {
      found_edges.at(2) = true;
    }
    if (edge == ed4) {
      found_edges.at(3) = true;
    }
    if (edge == ed5) {
      found_edges.at(4) = true;
    }
  }

  for (const bool& found : found_edges) {
    BOOST_CHECK(found);
  }
}

BOOST_AUTO_TEST_CASE(reduceGraph_algorithm_test2) {
  // Test that the reduced graph handles loops correctly
  //
  // 1 - 2
  // |   |
  // 4 - 3
  //
  // Should reduce this to
  //
  // 1 - 1
  //
  vector<Edge> edges{Edge(1, 2), Edge(2, 3), Edge(3, 4), Edge(1, 4)};

  unordered_map<votca::Index, GraphNode> nodes;
  for (votca::Index vertex = 1; vertex <= 4; ++vertex) {
    GraphNode gn;
    nodes[vertex] = gn;
  }

  Graph graph(edges, nodes);
  ReducedGraph reduced_graph = votca::tools::reduceGraph(graph);

  vector<Edge> edges2 = reduced_graph.getEdges();

  BOOST_CHECK_EQUAL(edges2.size(), 1);

  Edge ed_check(1, 1);
  BOOST_CHECK_EQUAL(edges2.at(0), ed_check);

  vector<vector<Edge>> edges3 = reduced_graph.expandEdge(ed_check);
  BOOST_CHECK_EQUAL(edges3.size(), 1);
  BOOST_CHECK_EQUAL(edges3.at(0).size(), 4);

  vector<bool> found_edge(4, false);

  for (const Edge& edge : edges3.at(0)) {
    if (edge == edges.at(0)) {
      found_edge.at(0) = true;
    }
    if (edge == edges.at(1)) {
      found_edge.at(1) = true;
    }
    if (edge == edges.at(2)) {
      found_edge.at(2) = true;
    }
    if (edge == edges.at(3)) {
      found_edge.at(3) = true;
    }
  }

  for (bool found : found_edge) {
    BOOST_CHECK(found);
  }
}

BOOST_AUTO_TEST_CASE(reduceGraph_algorithm_test3) {
  // Test that the reduced graph handles loop and edge correctly
  //
  // 1 - 2
  // |   |
  // 4 - 3 - 5 - 6
  //
  // Should reduce this to
  //
  // 3 - 3
  // 3 - 6
  //
  vector<Edge> edges{Edge(1, 2), Edge(2, 3), Edge(3, 4),
    Edge(1, 4), Edge(3, 5), Edge(5, 6)};

  unordered_map<votca::Index, GraphNode> nodes;
  for (votca::Index vertex = 1; vertex <= 6; ++vertex) {
    GraphNode gn;
    nodes[vertex] = gn;
  }

  Graph graph(edges, nodes);
  ReducedGraph reduced_graph = votca::tools::reduceGraph(graph);

  vector<Edge> edges2 = reduced_graph.getEdges();
  BOOST_CHECK_EQUAL(edges2.size(), 2);

  cout << "EDGE CHECK **************************************" << endl;

  Edge ed3_3(3, 3);
  Edge ed3_6(3, 6);

  bool found_edge3_3 = false;
  bool found_edge3_6 = false;
  for (Edge ed : edges2) {
    if (ed == ed3_3) {
      found_edge3_3 = true;
    } else if (ed == ed3_6) {
      found_edge3_6 = true;
    }
  }
  BOOST_CHECK(found_edge3_3);
  BOOST_CHECK(found_edge3_6);
}

BOOST_AUTO_TEST_CASE(reduceGraph_algorithm_test4) {
  // Test that the reduced graph handles loop and edge correctly
  //
  // 1 - 2       6 - 8
  // |   |       |   |
  // 4 - 3 - 5 - 7 - 9
  //
  // Should reduce this to
  //
  // 3 - 3
  // 3 - 7
  // 7 - 7
  //
  vector<Edge> edges{Edge(1, 2), Edge(2, 3), Edge(3, 4), Edge(1, 4),
    Edge(3, 5), Edge(5, 7), Edge(7, 6), Edge(6, 8),
    Edge(8, 9), Edge(9, 7)};

  unordered_map<votca::Index, GraphNode> nodes;
  for (votca::Index vertex = 1; vertex <= 9; ++vertex) {
    GraphNode gn;
    nodes[vertex] = gn;
  }

  Graph graph(edges, nodes);
  ReducedGraph reduced_graph = votca::tools::reduceGraph(graph);

  vector<votca::Index> junctions = reduced_graph.getJunctions();
  BOOST_CHECK_EQUAL(junctions.size(), 2);
  bool found_junction_3 = false;
  bool found_junction_7 = false;
  for (votca::Index junction : junctions) {
    if (junction == 3) {
      found_junction_3 = true;
    } else if (junction == 7) {
      found_junction_7 = true;
    }
  }
  BOOST_CHECK(found_junction_3);
  BOOST_CHECK(found_junction_7);

  vector<Edge> edges2 = reduced_graph.getEdges();
  BOOST_CHECK_EQUAL(edges2.size(), 3);

  Edge ed3_3(3, 3);
  Edge ed3_7(3, 7);
  Edge ed7_7(7, 7);

  bool found_edge3_3 = false;
  bool found_edge3_7 = false;
  bool found_edge7_7 = false;

  for (Edge ed : edges2) {
    if (ed == ed3_3) {
      found_edge3_3 = true;
    } else if (ed == ed3_7) {
      found_edge3_7 = true;
    } else if (ed == ed7_7) {
      found_edge7_7 = true;
    }
  }
  BOOST_CHECK(found_edge3_3);
  BOOST_CHECK(found_edge3_7);
  BOOST_CHECK(found_edge7_7);
}

BOOST_AUTO_TEST_CASE(reduceGraph_algorithm_test5) {

  // 4 - 5
  // |   |
  // 1 - 2
  //  \  |
  //    3
  Edge ed1(1, 2);
  Edge ed2(1, 4);
  Edge ed3(1, 3);
  Edge ed4(4, 5);
  Edge ed5(2, 5);
  Edge ed6(2, 3);

  vector<Edge> edges{ed1, ed2, ed3, ed4, ed5, ed6};

  unordered_map<votca::Index, GraphNode> nodes;
  for (votca::Index index = 0; index < 18; ++index) {
    GraphNode gn;
    nodes[index] = gn;
  }

  Graph g(edges, nodes);
  ReducedGraph reduced_g = votca::tools::reduceGraph(g);

  vector<votca::Index> junctions = reduced_g.getJunctions();
  BOOST_CHECK_EQUAL(junctions.size(), 2);
  bool found_junction_1 = false;
  bool found_junction_2 = false;
  for (votca::Index junction : junctions) {
    if (junction == 1) {
      found_junction_1 = true;
    } else if (junction == 2) {
      found_junction_2 = true;
    }
  }
  BOOST_CHECK(found_junction_1);
  BOOST_CHECK(found_junction_2);

  // Should end up with
  //  _ _
  // |   |
  // 1 - 2
  // |_ _|
  //

  vector<Edge> edges2 = reduced_g.getEdges();
  BOOST_CHECK_EQUAL(edges2.size(), 3);

  votca::Index edge_count1_2 = 0;
  for (const Edge& edge_temp : edges2) {
    if (edge_temp == ed1) {
      ++edge_count1_2;
    }
  }
  BOOST_CHECK_EQUAL(edge_count1_2, 3);
}

BOOST_AUTO_TEST_CASE(reduceGraph_algorithm_test6) {
  // Create edge
  //     2
  //     |
  // 0 - 1 - 3 - 9
  //     |       |
  //     4 - 5 - 6 -  7 - 8

  Edge ed1(0, 1);
  Edge ed2(1, 2);
  Edge ed3(1, 3);
  Edge ed4(1, 4);
  Edge ed5(4, 5);
  Edge ed6(5, 6);
  Edge ed7(6, 7);
  Edge ed8(7, 8);
  Edge ed9(6, 9);
  Edge ed10(3, 9);

  //
  // 10 - 11 - 12
  //
  // 13
  //
  // 14 - 15
  // |    |
  // 16 - 17
  //
  Edge ed11(10, 11);
  Edge ed12(11, 12);

  Edge ed13(14, 15);
  Edge ed14(15, 17);
  Edge ed15(16, 17);
  Edge ed16(14, 16);

  vector<Edge> edges{ed1, ed2,  ed3,  ed4,  ed5,  ed6,  ed7,  ed8,
    ed9, ed10, ed11, ed12, ed13, ed14, ed15, ed16};

  // Create Graph nodes
  unordered_map<votca::Index, GraphNode> nodes;
  for (votca::Index index = 0; index < 18; ++index) {
    GraphNode gn;
    nodes[index] = gn;
  }

  Graph g(edges, nodes);
  ReducedGraph reduced_g = votca::tools::reduceGraph(g);

  vector<votca::Index> junctions = reduced_g.getJunctions();
  BOOST_CHECK_EQUAL(junctions.size(), 2);

  bool found_junction_1 = false;
  bool found_junction_6 = false;
  for (votca::Index junction : junctions) {
    if (junction == 1) {
      found_junction_1 = true;
    } else if (junction == 6) {
      found_junction_6 = true;
    }
  }
  BOOST_CHECK(found_junction_1);
  BOOST_CHECK(found_junction_6);
  // Full Graph
  //
  //     2
  //     |
  // 0 - 1 - 3 - 9
  //     |       |
  //     4 - 5 - 6 -  7 - 8
  //
  // 10 - 11 - 12
  //
  // 13
  //
  // 14 - 15
  // |    |
  // 16 - 17
  //
  // Reduced Graph
  //
  //     2
  //     |
  // 0 - 1 -
  //     |   |
  //       - 6 - 8
  //
  // 10 - 12
  //
  // 13
  //
  // 14 -
  // |  |
  //  -
  //
  Edge ed0_1(0, 1);
  Edge ed1_2(1, 2);
  Edge ed1_6(1, 6);
  Edge ed6_8(6, 8);
  Edge ed10_12(10, 12);
  Edge ed14_14(14, 14);

  votca::Index edge_count0_1 = 0;
  votca::Index edge_count1_2 = 0;
  votca::Index edge_count1_6 = 0;
  votca::Index edge_count6_8 = 0;
  votca::Index edge_count10_12 = 0;
  votca::Index edge_count14_14 = 0;

  vector<Edge> edges2 = reduced_g.getEdges();

  BOOST_CHECK_EQUAL(edges2.size(), 7);

  for (Edge& edge : edges2) {
    if (edge == ed0_1) {
      ++edge_count0_1;
    }
    if (edge == ed1_2) {
      ++edge_count1_2;
    }
    if (edge == ed1_6) {
      ++edge_count1_6;
    }
    if (edge == ed6_8) {
      ++edge_count6_8;
    }
    if (edge == ed10_12) {
      ++edge_count10_12;
    }
    if (edge == ed14_14) {
      ++edge_count14_14;
    }
  }

  BOOST_CHECK_EQUAL(edge_count0_1, 1);
  BOOST_CHECK_EQUAL(edge_count1_2, 1);
  BOOST_CHECK_EQUAL(edge_count1_6, 2);
  BOOST_CHECK_EQUAL(edge_count6_8, 1);
  BOOST_CHECK_EQUAL(edge_count10_12, 1);
  BOOST_CHECK_EQUAL(edge_count14_14, 1);
}

BOOST_AUTO_TEST_CASE(explorebranch_test) {
  // Create edge
  //     2 - 10
  //     |
  // 0 - 1 - 3 - 9
  //     |       |
  //     4 - 5 - 6 -  7 - 8

  Edge ed1(0, 1);
  Edge ed2(1, 2);
  Edge ed3(1, 3);
  Edge ed4(1, 4);
  Edge ed5(4, 5);
  Edge ed6(5, 6);
  Edge ed7(6, 7);
  Edge ed8(7, 8);
  Edge ed9(6, 9);
  Edge ed10(3, 9);
  Edge ed11(2, 10);

  vector<Edge> edges{ed1, ed2, ed3, ed4, ed5, ed6, ed7, ed8, ed9, ed10, ed11};

  // Create Graph nodes
  unordered_map<votca::Index, GraphNode> nodes;
  for (votca::Index index = 0; index < 11; ++index) {
    GraphNode gn;
    nodes[index] = gn;
  }

  Graph g(edges, nodes);

  votca::Index starting_vertex = 1;
  set<Edge> branch_edges =
    votca::tools::exploreBranch(g, starting_vertex, ed3);

  // The following edges should be found in the branch
  //
  //     1 - 3 - 9
  //     |       |
  //     4 - 5 - 6 -  7 - 8
  //
  BOOST_CHECK_EQUAL(branch_edges.size(), 8);

  vector<bool> found_edges(8, false);
  votca::Index index = 0;
  for (const Edge& ed : branch_edges) {
    if (ed == ed3) {
      found_edges.at(index) = true;
    }
    if (ed == ed4) {
      found_edges.at(index) = true;
    }
    if (ed == ed5) {
      found_edges.at(index) = true;
    }
    if (ed == ed6) {
      found_edges.at(index) = true;
    }
    if (ed == ed7) {
      found_edges.at(index) = true;
    }
    if (ed == ed8) {
      found_edges.at(index) = true;
    }
    if (ed == ed9) {
      found_edges.at(index) = true;
    }
    if (ed == ed10) {
      found_edges.at(index) = true;
    }
    ++index;
  }

  for (const bool& found : found_edges) {
    BOOST_CHECK(found);
  }
}

BOOST_AUTO_TEST_CASE(structureid_test) {

  // Create edge
  Edge ed(0, 1);
  Edge ed1(1, 2);
  Edge ed2(2, 3);
  Edge ed3(3, 4);
  Edge ed4(4, 5);
  Edge ed5(5, 0);
  Edge ed6(3, 6);

  vector<Edge> edges;
  edges.push_back(ed);
  edges.push_back(ed1);
  edges.push_back(ed2);
  edges.push_back(ed3);
  edges.push_back(ed4);
  edges.push_back(ed5);
  edges.push_back(ed6);

  // Create Graph nodes
  GraphNode gn1;
  GraphNode gn2;
  GraphNode gn3;
  GraphNode gn4;
  GraphNode gn5;
  GraphNode gn6;
  GraphNode gn7;

  unordered_map<votca::Index, GraphNode> nodes;
  nodes[0] = gn1;
  nodes[1] = gn2;
  nodes[2] = gn3;
  nodes[3] = gn4;
  nodes[4] = gn5;
  nodes[5] = gn6;
  nodes[6] = gn7;

  Graph g(edges, nodes);

  string structId = votca::tools::findStructureId<GraphDistVisitor>(g).get();

  string answer = "Dist=0;Dist=1;Dist=1;Dist=1;Dist=2;Dist=2;Dist=3;";
  BOOST_CHECK_EQUAL(structId, answer);
}
/*
   BOOST_AUTO_TEST_CASE(find_canonized_sequence_test) {
   {

// Graph A       Graph B
// 1a - 2a - 3a     4a - 3a - 2a
//       |   |            |   |
//      5b - 4b          1b - 5b
///
// Create edge
Edge ed1_A(1, 2);
Edge ed2_A(2, 3);
Edge ed3_A(3, 4);
Edge ed4_A(4, 5);
Edge ed5_A(5, 2);

vector<Edge> edges_A;
edges_A.push_back(ed1_A);
edges_A.push_back(ed2_A);
edges_A.push_back(ed3_A);
edges_A.push_back(ed4_A);
edges_A.push_back(ed5_A);

// Create Graph nodes
GraphNode gn1_A;
GraphNode gn2_A;
GraphNode gn3_A;
GraphNode gn4_A;
GraphNode gn5_A;

gn1_A.addStr("label", "a");
gn2_A.addStr("label", "a");
gn3_A.addStr("label", "a");
gn4_A.addStr("label", "b");
gn5_A.addStr("label", "b");

unordered_map<votca::Index, GraphNode> nodes_A;
nodes_A[1] = gn1_A;
nodes_A[2] = gn2_A;
nodes_A[3] = gn3_A;
nodes_A[4] = gn4_A;
nodes_A[5] = gn5_A;

Graph gA(edges_A, nodes_A);
std::vector<votca::Index> sequence_A;
std::string structId_A =
votca::tools::findCanonizedSequence(gA, sequence_A);

// Graph A       Graph B
// 1a - 2a - 3a     4a - 3a - 2a
//       |   |            |   |
//      5b - 4b          1b - 5b
///
// Create edge
Edge ed1_B(4, 3);
Edge ed2_B(3, 2);
Edge ed3_B(2, 5);
Edge ed4_B(5, 1);
Edge ed5_B(1, 3);

vector<Edge> edges_B;
edges_B.push_back(ed1_B);
edges_B.push_back(ed2_B);
edges_B.push_back(ed3_B);
edges_B.push_back(ed4_B);
edges_B.push_back(ed5_B);

// Create Graph nodes
GraphNode gn1_B;
GraphNode gn2_B;
GraphNode gn3_B;
GraphNode gn4_B;
GraphNode gn5_B;

gn1_B.addStr("label", "b");
gn2_B.addStr("label", "a");
gn3_B.addStr("label", "a");
gn4_B.addStr("label", "a");
gn5_B.addStr("label", "b");

unordered_map<votca::Index, GraphNode> nodes_B;
nodes_B[1] = gn1_B;
nodes_B[2] = gn2_B;
nodes_B[3] = gn3_B;
nodes_B[4] = gn4_B;
nodes_B[5] = gn5_B;

Graph gB(edges_B, nodes_B);
vector<votca::Index> sequence_B;
string structId_B = votca::tools::findCanonizedSequence(gB, sequence_B);

cout << "struct ID A " << structId_A << endl;
cout << "struct ID B " << structId_B << endl;

cout << "Print Sequence" << endl;
cout << "A     B    Cannonized A     Cannonzied B     Degree A   Degree B"
<< endl;
cout << "================================================================"
<< endl;
for (int i = 0; i < 5; ++i) {
  cout << sequence_A.at(i) << "     " << sequence_B.at(i);
  cout << "          ";
  cout << nodes_A.at(sequence_A.at(i)).getStr("label");
  cout << "               ";
  cout << nodes_B.at(sequence_B.at(i)).getStr("label");
  cout << "              ";
  cout << gA.getDegree(sequence_A.at(i));
  cout << "          ";
  cout << gB.getDegree(sequence_B.at(i));
  cout << endl;
}
string answer = "Dist0labelaDist1labelaDist1labelaDist1labelbDist2labelb";
BOOST_CHECK_EQUAL(structId_A, answer);
BOOST_CHECK_EQUAL(structId_B, answer);

BOOST_CHECK_EQUAL(sequence_A.at(0), 2);
BOOST_CHECK_EQUAL(sequence_A.at(1), 5);
BOOST_CHECK_EQUAL(sequence_A.at(2), 1);
BOOST_CHECK_EQUAL(sequence_A.at(3), 3);
BOOST_CHECK_EQUAL(sequence_A.at(4), 4);

BOOST_CHECK_EQUAL(sequence_B.at(0), 3);
BOOST_CHECK_EQUAL(sequence_B.at(1), 1);
BOOST_CHECK_EQUAL(sequence_B.at(2), 4);
BOOST_CHECK_EQUAL(sequence_B.at(3), 2);
BOOST_CHECK_EQUAL(sequence_B.at(4), 5);

BOOST_CHECK_EQUAL(gA.getDegree(sequence_A.at(0)), 3);
BOOST_CHECK_EQUAL(gA.getDegree(sequence_A.at(1)), 2);
BOOST_CHECK_EQUAL(gA.getDegree(sequence_A.at(2)), 1);
BOOST_CHECK_EQUAL(gA.getDegree(sequence_A.at(3)), 2);
BOOST_CHECK_EQUAL(gA.getDegree(sequence_A.at(4)), 2);

BOOST_CHECK_EQUAL(gB.getDegree(sequence_B.at(0)), 3);
BOOST_CHECK_EQUAL(gB.getDegree(sequence_B.at(1)), 2);
BOOST_CHECK_EQUAL(gB.getDegree(sequence_B.at(2)), 1);
BOOST_CHECK_EQUAL(gB.getDegree(sequence_B.at(3)), 2);
BOOST_CHECK_EQUAL(gB.getDegree(sequence_B.at(4)), 2);
}
}*/
BOOST_AUTO_TEST_SUITE_END()
