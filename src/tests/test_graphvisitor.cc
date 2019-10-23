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

#define BOOST_TEST_MODULE graphvisitor_test
#include <boost/test/unit_test.hpp>
#include <unordered_map>
#include <vector>
#include <votca/tools/graph.h>
#include <votca/tools/graphnode.h>
#include <votca/tools/graphvisitor.h>

using namespace std;
using namespace votca::tools;

class GraphVisitorTest : public GraphVisitor {
 private:
  void addEdges_(const Graph&, long int) override {
    throw runtime_error("Undefined method.");
  }
  Edge getEdge_() override {
    throw runtime_error("Undefined method.");
  }
};

BOOST_AUTO_TEST_SUITE(graphvisitor_test)

BOOST_AUTO_TEST_CASE(constructor_test) { GraphVisitorTest gv; }

BOOST_AUTO_TEST_CASE(basic_test) {
  // Create edge
  Edge ed(0, 1);
  vector<Edge> edges;
  edges.push_back(ed);

  // Create Graph nodes
  GraphNode gn1;
  GraphNode gn2;

  unordered_map<long int, GraphNode> nodes;
  nodes[0] = gn1;
  nodes[1] = gn2;

  Graph g(edges, nodes);

  GraphVisitorTest gv;

  BOOST_CHECK(gv.queEmpty());

  BOOST_CHECK_THROW(gv.exec(g, ed), runtime_error);
  BOOST_CHECK_THROW(gv.initialize(g), runtime_error);
  // No exception should be thrown at this point
  gv.exec(g, ed);
  // Error because no nextEdge function ptr passed in
  BOOST_CHECK_THROW(gv.nextEdge(g), runtime_error);

  BOOST_CHECK_EQUAL(gv.getStartingVertex(), 0);
  gv.setStartingVertex(2);
  BOOST_CHECK_EQUAL(gv.getStartingVertex(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
