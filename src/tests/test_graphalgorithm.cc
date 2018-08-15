/*
 *            Copyright 2009-2018 The VOTCA Development Team
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
#include <unordered_map>
#include <vector>
#include <votca/tools/graph.h>
#include <votca/tools/graph_bf_visitor.h>
#include <votca/tools/graphalgorithm.h>
#include <votca/tools/graphdistvisitor.h>
#include <votca/tools/graphnode.h>

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
    BOOST_CHECK_THROW(gb_v.exec(g, ed), runtime_error);

    bool single_n = singleNetwork(g, gb_v);
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

BOOST_AUTO_TEST_CASE(decoupleIsolatedSubGraphs_algorithm_test) {
  {
    // In this test we add two nodes and an edge describing
    // their connection

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

    auto subGraphs = decoupleIsolatedSubGraphs(g);

    BOOST_CHECK_EQUAL(subGraphs.size(), 1);
    BOOST_CHECK_EQUAL(*(subGraphs.at(0)), g);
  }

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
    vector<std::shared_ptr<Graph>> sub_graphs = decoupleIsolatedSubGraphs(g);
    BOOST_CHECK_EQUAL(sub_graphs.size(), 2);

    // Create sub graphs to compare with the graphs stored in 'sub_graphs'

    vector<Edge> sub_graph1_edges;
    sub_graph1_edges.push_back(ed1);
    sub_graph1_edges.push_back(ed2);
    sub_graph1_edges.push_back(ed3);
    sub_graph1_edges.push_back(ed4);

    vector<Edge> sub_graph2_edges;
    sub_graph2_edges.push_back(ed5);

    unordered_map<int, GraphNode> sub_graph1_nodes;
    sub_graph1_nodes[0] = gn1;
    sub_graph1_nodes[1] = gn2;
    sub_graph1_nodes[2] = gn3;
    sub_graph1_nodes[3] = gn4;
    sub_graph1_nodes[4] = gn5;

    unordered_map<int, GraphNode> sub_graph2_nodes;
    sub_graph2_nodes[0] = gn6;
    sub_graph2_nodes[1] = gn7;

    Graph sub_graph1(sub_graph1_edges, sub_graph1_nodes);
    Graph sub_graph2(sub_graph2_edges, sub_graph2_nodes);

    // Cycle sub_graphs
    bool sub_graph1_found = false;
    bool sub_graph2_found = false;

    for (auto graph_sub_it = sub_graphs.begin();
         graph_sub_it != sub_graphs.end(); ++graph_sub_it) {
      if (**graph_sub_it == sub_graph1) sub_graph1_found = true;
      if (**graph_sub_it == sub_graph2) sub_graph2_found = true;
    }
    BOOST_CHECK(sub_graph1_found);
    BOOST_CHECK(sub_graph2_found);
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

    auto structId = findStructureId<GraphDistVisitor>(g);

    string answer = "Dist0Dist1Dist1Dist1Dist2Dist2Dist3";
    BOOST_CHECK_EQUAL(structId, answer);
  }
}

BOOST_AUTO_TEST_SUITE_END()
