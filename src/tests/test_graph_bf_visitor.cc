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

#define BOOST_TEST_MODULE graph_bf_visitor_test
#include <boost/test/unit_test.hpp>
#include <unordered_map>
#include <vector>
#include <votca/tools/graph.h>
#include <votca/tools/graph_bf_visitor.h>
#include <votca/tools/graphnode.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(graph_bf_visitor_test)

BOOST_AUTO_TEST_CASE(constructor_test) { Graph_BF_Visitor gb_v; }

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

  Graph_BF_Visitor gb_v;
  BOOST_CHECK(gb_v.queEmpty());
  BOOST_CHECK_THROW(gb_v.exec(g, ed), runtime_error);
  // Default starts with node index 0
  gb_v.initialize(g);
  BOOST_CHECK_EQUAL(gb_v.queEmpty(), false);
  // No exception should be thrown at this point
  Edge ed1 = gb_v.nextEdge(g);
  BOOST_CHECK_EQUAL(ed, ed1);
  gb_v.exec(g, ed1);
  BOOST_CHECK(gb_v.queEmpty());

  // Show which vertices have been explored
  auto exploredV = gb_v.getExploredVertices();

  bool v0 = false;
  bool v1 = false;
  for (auto ver : exploredV) {
    if (ver == 0) v0 = true;
    if (ver == 1) v1 = true;
  }

  // Both vertices should have been explored
  BOOST_CHECK(v0);
  BOOST_CHECK(v1);
}

BOOST_AUTO_TEST_CASE(basic_test2) {

  // This tests demonstrates that this particular graph visitor
  // explores the graph in a breadth first, and first in first out order
  //
  // 0 -> 1 -> 2
  // |
  // v
  // 3
  // |
  // v
  // 4
  //

  // Create edge
  Edge ed(0, 1);
  Edge ed1(1, 2);
  Edge ed2(0, 3);
  Edge ed3(3, 4);

  vector<Edge> edges;
  edges.push_back(ed);
  edges.push_back(ed1);
  edges.push_back(ed2);
  edges.push_back(ed3);

  // Create Graph nodes
  GraphNode gn1;
  GraphNode gn2;
  GraphNode gn3;
  GraphNode gn4;
  GraphNode gn5;

  unordered_map<int, GraphNode> nodes;
  nodes[0] = gn1;
  nodes[1] = gn2;
  nodes[2] = gn3;
  nodes[3] = gn4;
  nodes[4] = gn5;

  Graph g(edges, nodes);

  Graph_BF_Visitor gb_v;
  BOOST_CHECK(gb_v.queEmpty());
  BOOST_CHECK_THROW(gb_v.exec(g, ed), runtime_error);
  // Default starts with node index 0
  gb_v.initialize(g);
  BOOST_CHECK_EQUAL(gb_v.queEmpty(), false);
  // No exception should be thrown at this point

  // First two edges found should be ed and ed2
  vector<Edge> temp;
  Edge ed5 = gb_v.nextEdge(g);
  temp.push_back(ed5);
  gb_v.exec(g, ed5);
  ed5 = gb_v.nextEdge(g);
  temp.push_back(ed5);

  bool found_ed = false;
  bool found_ed2 = false;
  for (auto temp_ed : temp) {
    if (temp_ed == ed) found_ed = true;
    if (temp_ed == ed2) found_ed2 = true;
  }
  BOOST_CHECK(found_ed);
  BOOST_CHECK(found_ed2);

  // Next two edges are ed1 and ed3 they are equal distance apart from the
  // starting vertex so the order does not matter.
  temp.clear();

  gb_v.exec(g, ed5);
  ed5 = gb_v.nextEdge(g);
  temp.push_back(ed5);
  gb_v.exec(g, ed5);
  ed5 = gb_v.nextEdge(g);
  temp.push_back(ed5);

  bool found_ed1 = false;
  bool found_ed3 = false;
  for (auto temp_ed : temp) {
    if (temp_ed == ed1) found_ed1 = true;
    if (temp_ed == ed3) found_ed3 = true;
  }
  BOOST_CHECK(found_ed1);
  BOOST_CHECK(found_ed3);

  gb_v.exec(g, ed5);
  BOOST_CHECK(gb_v.queEmpty());

  // Show which vertices have been explored
  auto exploredV = gb_v.getExploredVertices();

  bool v0 = false;
  bool v1 = false;
  bool v2 = false;
  bool v3 = false;
  bool v4 = false;
  for (auto ver : exploredV) {
    if (ver == 0) v0 = true;
    if (ver == 1) v1 = true;
    if (ver == 2) v2 = true;
    if (ver == 3) v3 = true;
    if (ver == 4) v4 = true;
  }

  // All 5 vertices should have been explored
  BOOST_CHECK(v0);
  BOOST_CHECK(v1);
  BOOST_CHECK(v2);
  BOOST_CHECK(v3);
  BOOST_CHECK(v4);
}

BOOST_AUTO_TEST_SUITE_END()
