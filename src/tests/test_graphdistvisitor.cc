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

#define BOOST_TEST_MODULE graphdistvisitor_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <votca/tools/graph.h>
#include <votca/tools/graphdistvisitor.h>
#include <votca/tools/graphnode.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(graphdistvisitor_test)

BOOST_AUTO_TEST_CASE(constructor_test) { GraphDistVisitor gdv; }

BOOST_AUTO_TEST_CASE(basic_test) {
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

  GraphDistVisitor gdv;
  BOOST_CHECK(gdv.queEmpty());

  BOOST_CHECK_THROW(gdv.exec(g, ed), runtime_error);
  // Default starts with node index 0
  gdv.initialize(g);
  BOOST_CHECK_EQUAL(gdv.queEmpty(), false);
  // No exception should be thrown at this point
  Edge ed1 = gdv.nextEdge(g);
  BOOST_CHECK_EQUAL(ed, ed1);
  gdv.exec(g, ed1);
  BOOST_CHECK(gdv.queEmpty());
  GraphNode gn3 = g.getNode(0);
  int dist = gn3.getInt("Dist");
  BOOST_CHECK_EQUAL(dist, 0);
  GraphNode gn4 = g.getNode(1);
  dist = gn4.getInt("Dist");
  BOOST_CHECK_EQUAL(dist, 1);
}

BOOST_AUTO_TEST_CASE(basic_test2) {
  // Create edges
  //
  // 0 - 1 - 2
  // |       |
  // 6 - 4 - 3
  //     |
  //     5
  //
  Edge ed(0, 1);
  Edge ed1(1, 2);
  Edge ed2(2, 3);
  Edge ed3(3, 4);
  Edge ed4(4, 5);
  Edge ed5(4, 6);
  Edge ed6(6, 0);

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

  GraphDistVisitor gdv;
  BOOST_CHECK(gdv.queEmpty());

  BOOST_CHECK_THROW(gdv.exec(g, ed), runtime_error);
  // Default starts with node index 0
  gdv.initialize(g);
  BOOST_CHECK_EQUAL(gdv.queEmpty(), false);
  // No exception should be thrown at this point

  // First two edges that should be explored are edges ed and ed6
  vector<Edge> temp;
  Edge ed7 = gdv.nextEdge(g);
  temp.push_back(ed7);
  gdv.exec(g, ed7);
  ed7 = gdv.nextEdge(g);
  gdv.exec(g, ed7);
  temp.push_back(ed7);

  bool found_ed = false;
  bool found_ed6 = false;
  for (auto ed_temp : temp) {
    if (ed_temp == ed) found_ed = true;
    if (ed_temp == ed6) found_ed6 = true;
  }
  BOOST_CHECK(found_ed);
  BOOST_CHECK(found_ed6);

  // Explore the whole graph
  while (!gdv.queEmpty()) {
    ed7 = gdv.nextEdge(g);
    gdv.exec(g, ed7);
  }

  // Now let's check that the distances for each node
  // are correct. The distances are counted from the starting
  // node which by default is node 0.

  GraphNode gn8 = g.getNode(0);
  int dist = gn8.getInt("Dist");
  BOOST_CHECK_EQUAL(dist, 0);

  gn8 = g.getNode(1);
  dist = gn8.getInt("Dist");
  BOOST_CHECK_EQUAL(dist, 1);

  gn8 = g.getNode(2);
  dist = gn8.getInt("Dist");
  BOOST_CHECK_EQUAL(dist, 2);

  gn8 = g.getNode(3);
  dist = gn8.getInt("Dist");
  BOOST_CHECK_EQUAL(dist, 3);

  gn8 = g.getNode(4);
  dist = gn8.getInt("Dist");
  BOOST_CHECK_EQUAL(dist, 2);

  gn8 = g.getNode(5);
  dist = gn8.getInt("Dist");
  BOOST_CHECK_EQUAL(dist, 3);

  gn8 = g.getNode(6);
  dist = gn8.getInt("Dist");
  BOOST_CHECK_EQUAL(dist, 1);
}

BOOST_AUTO_TEST_SUITE_END()
