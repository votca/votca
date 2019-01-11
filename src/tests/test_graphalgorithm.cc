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
#include <memory>
#include <unordered_map>
#include <vector>
#include <votca/tools/graph.h>
#include <votca/tools/reducedgraph.h>
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
    //gb_v.exec(g,ed);
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

    Graph g(edges,nodes);
    ReducedGraph reduced_g = reduceGraph(g);

    auto edges2 = reduced_g.getEdges();
    BOOST_CHECK_EQUAL(edges2.size(),5);
    vector<bool> found_edges(5,false);
    for(auto edge : edges2 ){
      if(edge==ed1) found_edges.at(0) = true;
      if(edge==ed2) found_edges.at(1) = true;
      if(edge==ed3) found_edges.at(2) = true;
      if(edge==ed4) found_edges.at(3) = true;
      if(edge==ed5) found_edges.at(4) = true;
    }

    for(auto found : found_edges){
      BOOST_CHECK(found);
    }
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

    Edge ed13(14,15);
    Edge ed14(15,16);
    Edge ed15(15,17);
    Edge ed16(14,16);

    vector<Edge> edges{ ed1, ed2, ed3, ed4, ed5, ed6, ed7, ed8, ed9, ed10,
      ed11, ed12, ed13, ed14, ed15, ed16};

    // Create Graph nodes
    unordered_map<int, GraphNode> nodes;
    for(int index = 0; index<18;++index){
      GraphNode gn;
      nodes[index] = gn;
    }

    Graph g(edges,nodes);
    ReducedGraph reduced_g = reduceGraph(g);

    auto edges2 = reduced_g.getEdges();
    BOOST_CHECK_EQUAL(edges2.size(),5);
    vector<bool> found_edges(5,false);
    for(auto edge : edges2 ){
      if(edge==ed1) found_edges.at(0) = true;
      if(edge==ed2) found_edges.at(1) = true;
      if(edge==ed3) found_edges.at(2) = true;
      if(edge==ed4) found_edges.at(3) = true;
      if(edge==ed5) found_edges.at(4) = true;
    }

    for(auto found : found_edges){
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

    auto structId = findStructureId<GraphDistVisitor>(g);

    string answer = "Dist0Dist1Dist1Dist1Dist2Dist2Dist3";
    BOOST_CHECK_EQUAL(structId, answer);
  }
}

BOOST_AUTO_TEST_SUITE_END()
