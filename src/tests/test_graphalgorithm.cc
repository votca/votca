/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <boost/test/unit_test.hpp>
#include <memory>
#include <unordered_map>
#include <vector>
#include <votca/tools/graph.h>
#include <votca/tools/graph_bf_visitor.h>
#include <votca/tools/graphalgorithm.h>
#include <votca/tools/graphdistvisitor.h>
#include <votca/tools/graphnode.h>
#include <votca/tools/reducedgraph.h>

using namespace std;
using namespace votca::tools;

using namespace boost;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(graphalgorithm_test)

BOOST_AUTO_TEST_CASE(single_network_algorithm_test) {
  {
    // In this test we add two nodes and an edge describing
    // their connection thus the singleNetwork function will
    // return true.

    // Create edge
    Edge ed(0, 1);
    vector<Edge> edges;
    edges.push_back(ed);

    // Create Graph nodes
    GraphNode gn1;
    GraphNode gn2;

    unordered_map<int, GraphNode> nodes;
    nodes[0] = gn1;
    nodes[1] = gn2;

    Graph g(edges, nodes);

    Graph_BF_Visitor gb_v;

    BOOST_CHECK(gb_v.queEmpty());
    // gb_v.exec(g,ed);
    BOOST_CHECK_THROW(gb_v.exec(g, ed), runtime_error);

    bool single_n = singleNetwork(g, gb_v);
    cerr << "is single network " << single_n << endl;
    BOOST_CHECK(single_n);
    BOOST_CHECK(gb_v.queEmpty());

    Graph_BF_Visitor gb_v2;
    gb_v2.setStartingVertex(2);
    BOOST_CHECK_THROW(singleNetwork(g, gb_v2), invalid_argument);
  }

  {

    // In this test we add 3 nodes but only one edge
    // this means that one of the nodes will not be
    // attached to the other two. Thus the singleNetwork
    // function should return false.

    // Create edge
    Edge ed(0, 1);
    vector<Edge> edges;
    edges.push_back(ed);

    // Create Graph nodes
    GraphNode gn1;
    GraphNode gn2;
    GraphNode gn3;

    unordered_map<int, GraphNode> nodes;
    nodes[0] = gn1;
    nodes[1] = gn2;
    nodes[2] = gn3;

    Graph g(edges, nodes);

    Graph_BF_Visitor gb_v;

    BOOST_CHECK(gb_v.queEmpty());
    BOOST_CHECK_THROW(gb_v.exec(g, ed), runtime_error);

    bool single_n = singleNetwork(g, gb_v);
    BOOST_CHECK(!single_n);
    BOOST_CHECK(gb_v.queEmpty());
  }
}

BOOST_AUTO_TEST_CASE(reduceGraph_algorithm_test) {

  {
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

    unordered_map<int, GraphNode> nodes;
    nodes[0] = gn1;
    nodes[1] = gn2;
    nodes[2] = gn3;
    nodes[3] = gn4;
    nodes[4] = gn5;
    nodes[5] = gn6;
    nodes[6] = gn7;

    Graph g(edges, nodes);
    ReducedGraph reduced_g = reduceGraph(g);

    vector<Edge> edges2 = reduced_g.getEdges();
    BOOST_CHECK_EQUAL(edges2.size(), 5);
    vector<bool> found_edges(5, false);
    for (Edge& edge : edges2) {
      if (edge == ed1) found_edges.at(0) = true;
      if (edge == ed2) found_edges.at(1) = true;
      if (edge == ed3) found_edges.at(2) = true;
      if (edge == ed4) found_edges.at(3) = true;
      if (edge == ed5) found_edges.at(4) = true;
    }

    for (const bool& found : found_edges) {
      BOOST_CHECK(found);
    }
  }

  {
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

    unordered_map<int, GraphNode> nodes;
    for (int vertex = 1; vertex <= 4; ++vertex) {
      GraphNode gn;
      nodes[vertex] = gn;
    }

    Graph graph(edges, nodes);
    ReducedGraph reduced_graph = reduceGraph(graph);

    vector<Edge> edges2 = reduced_graph.getEdges();

    BOOST_CHECK_EQUAL(edges2.size(), 1);

    Edge ed_check(1, 1);
    BOOST_CHECK_EQUAL(edges2.at(0), ed_check);

    vector<vector<Edge>> edges3 = reduced_graph.expandEdge(ed_check);
    BOOST_CHECK_EQUAL(edges3.size(), 1);
    BOOST_CHECK_EQUAL(edges3.at(0).size(), 4);

    vector<bool> found_edge(4, false);

    for (const Edge& edge : edges3.at(0)) {
      if (edge == edges.at(0)) found_edge.at(0) = true;
      if (edge == edges.at(1)) found_edge.at(1) = true;
      if (edge == edges.at(2)) found_edge.at(2) = true;
      if (edge == edges.at(3)) found_edge.at(3) = true;
    }

    for (bool found : found_edge) {
      BOOST_CHECK(found);
    }
  }

  {
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

    unordered_map<int, GraphNode> nodes;
    for (int vertex = 1; vertex <= 6; ++vertex) {
      GraphNode gn;
      nodes[vertex] = gn;
    }

    Graph graph(edges, nodes);
    ReducedGraph reduced_graph = reduceGraph(graph);

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

  {
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

    unordered_map<int, GraphNode> nodes;
    for (int vertex = 1; vertex <= 9; ++vertex) {
      GraphNode gn;
      nodes[vertex] = gn;
    }

    Graph graph(edges, nodes);
    ReducedGraph reduced_graph = reduceGraph(graph);

    vector<int> junctions = reduced_graph.getJunctions();
    BOOST_CHECK_EQUAL(junctions.size(), 2);
    bool found_junction_3 = false;
    bool found_junction_7 = false;
    for (int junction : junctions) {
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

  {

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

    unordered_map<int, GraphNode> nodes;
    for (int index = 0; index < 18; ++index) {
      GraphNode gn;
      nodes[index] = gn;
    }

    Graph g(edges, nodes);
    ReducedGraph reduced_g = reduceGraph(g);

    vector<int> junctions = reduced_g.getJunctions();
    BOOST_CHECK_EQUAL(junctions.size(), 2);
    bool found_junction_1 = false;
    bool found_junction_2 = false;
    for (int junction : junctions) {
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

    int edge_count1_2 = 0;
    for (const Edge& edge_temp : edges2) {
      if (edge_temp == ed1) ++edge_count1_2;
    }
    BOOST_CHECK_EQUAL(edge_count1_2, 3);
  }

  {
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
    unordered_map<int, GraphNode> nodes;
    for (int index = 0; index < 18; ++index) {
      GraphNode gn;
      nodes[index] = gn;
    }

    Graph g(edges, nodes);
    ReducedGraph reduced_g = reduceGraph(g);

    vector<int> junctions = reduced_g.getJunctions();
    BOOST_CHECK_EQUAL(junctions.size(), 2);

    bool found_junction_1 = false;
    bool found_junction_6 = false;
    for (int junction : junctions) {
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

    int edge_count0_1 = 0;
    int edge_count1_2 = 0;
    int edge_count1_6 = 0;
    int edge_count6_8 = 0;
    int edge_count10_12 = 0;
    int edge_count14_14 = 0;

    vector<Edge> edges2 = reduced_g.getEdges();

    BOOST_CHECK_EQUAL(edges2.size(), 7);

    for (Edge& edge : edges2) {
      if (edge == ed0_1) ++edge_count0_1;
      if (edge == ed1_2) ++edge_count1_2;
      if (edge == ed1_6) ++edge_count1_6;
      if (edge == ed6_8) ++edge_count6_8;
      if (edge == ed10_12) ++edge_count10_12;
      if (edge == ed14_14) ++edge_count14_14;
    }

    BOOST_CHECK_EQUAL(edge_count0_1, 1);
    BOOST_CHECK_EQUAL(edge_count1_2, 1);
    BOOST_CHECK_EQUAL(edge_count1_6, 2);
    BOOST_CHECK_EQUAL(edge_count6_8, 1);
    BOOST_CHECK_EQUAL(edge_count10_12, 1);
    BOOST_CHECK_EQUAL(edge_count14_14, 1);
  }
}

BOOST_AUTO_TEST_CASE(explorebranch_test) {
  {
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
    unordered_map<int, GraphNode> nodes;
    for (int index = 0; index < 11; ++index) {
      GraphNode gn;
      nodes[index] = gn;
    }

    Graph g(edges, nodes);

    int starting_vertex = 1;
    set<Edge> branch_edges = exploreBranch(g, starting_vertex, ed3);

    // The following edges should be found in the branch
    //
    //     1 - 3 - 9
    //     |       |
    //     4 - 5 - 6 -  7 - 8
    //
    BOOST_CHECK_EQUAL(branch_edges.size(), 8);

    vector<bool> found_edges(8, false);
    int index = 0;
    for (const Edge& ed : branch_edges) {
      if (ed == ed3) found_edges.at(index) = true;
      if (ed == ed4) found_edges.at(index) = true;
      if (ed == ed5) found_edges.at(index) = true;
      if (ed == ed6) found_edges.at(index) = true;
      if (ed == ed7) found_edges.at(index) = true;
      if (ed == ed8) found_edges.at(index) = true;
      if (ed == ed9) found_edges.at(index) = true;
      if (ed == ed10) found_edges.at(index) = true;
      ++index;
    }

    for (const bool& found : found_edges) {
      BOOST_CHECK(found);
    }
  }
}

BOOST_AUTO_TEST_CASE(structureid_test) {
  {

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

    unordered_map<int, GraphNode> nodes;
    nodes[0] = gn1;
    nodes[1] = gn2;
    nodes[2] = gn3;
    nodes[3] = gn4;
    nodes[4] = gn5;
    nodes[5] = gn6;
    nodes[6] = gn7;

    Graph g(edges, nodes);

    string structId = findStructureId<GraphDistVisitor>(g);

    string answer = "Dist0Dist1Dist1Dist1Dist2Dist2Dist3";
    BOOST_CHECK_EQUAL(structId, answer);
  }
}

BOOST_AUTO_TEST_SUITE_END()
