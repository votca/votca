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

#define BOOST_TEST_MODULE graphbasicvisitor_test
#include <boost/test/unit_test.hpp>
#include <vector>
#include <unordered_map>
#include <votca/tools/graph.h>
#include <votca/tools/graphnode.h>
#include <votca/tools/graphbasicvisitor.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(graphbasicvisitor_test)

BOOST_AUTO_TEST_CASE(constructor_test){
  GraphBasicVisitor gb_v;
}

BOOST_AUTO_TEST_CASE(basic_test){
  // Create edge
  Edge ed(0,1);
  vector<Edge> edges;
  edges.push_back(ed);
  
  // Create Graph nodes
  GraphNode gn1;
  GraphNode gn2;
  
  unordered_map<int,GraphNode> nodes;
  nodes[0] = gn1;
  nodes[1] = gn2;

  Graph g(edges,nodes);

  GraphBasicVisitor gb_v;
  BOOST_CHECK(gb_v.queEmpty());
  BOOST_CHECK_THROW(gb_v.exec(g,ed),runtime_error);
  // Default starts with node index 0
  gb_v.startingVertex(g);
  BOOST_CHECK_EQUAL(gb_v.queEmpty(),false);
  // No exception should be thrown at this point
  Edge ed1 = gb_v.nextEdge(g);
  BOOST_CHECK_EQUAL(ed,ed1);
  gb_v.exec(g,ed1);
  BOOST_CHECK(gb_v.queEmpty());
}

BOOST_AUTO_TEST_SUITE_END()
