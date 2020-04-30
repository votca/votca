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

#define BOOST_TEST_MODULE branch_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "../libtools/branch.h"
#include "votca/tools/graphnode.h"
#include "votca/tools/reducededge.h"
#include "votca/tools/reducedgraph.h"

using namespace std;
using namespace votca;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(branch_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  vector<ReducedEdge> vec_ed;
  ReducedEdge edge(std::vector<Index>{0, 1, 2});
  vec_ed.push_back(edge);

  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;

  ReducedGraph graph(vec_ed, m_gn);

  Index starting_vertex = 0;
  Branch branch(edge, starting_vertex, graph);
}

BOOST_AUTO_TEST_CASE(terminal_and_source) {
  vector<ReducedEdge> vec_ed;
  ReducedEdge edge(std::vector<Index>{0, 1, 2});
  vec_ed.push_back(edge);

  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;

  ReducedGraph graph(vec_ed, m_gn);

  Index starting_vertex = 0;
  Branch branch(edge, starting_vertex, graph);

  BOOST_CHECK(branch.getSource() == starting_vertex);
  Index ending_vertex = 2;
  BOOST_CHECK(branch.getTerminal() == ending_vertex);
}

BOOST_AUTO_TEST_CASE(string_id) {
  vector<ReducedEdge> vec_ed;
  ReducedEdge edge(std::vector<Index>{0, 1, 2});
  vec_ed.push_back(edge);

  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  gn.add(string("name"), string("a"));
  gn1.add(string("name"), string("c"));
  gn2.add(string("name"), string("b"));

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;

  ReducedGraph graph(vec_ed, m_gn);

  Index starting_vertex = 0;
  Branch branch(edge, starting_vertex, graph);

  BOOST_TEST(branch.getBranchStringId() == string("name=a;name=c;name=b;"));
}

BOOST_AUTO_TEST_CASE(reverse_sequence) {
  vector<ReducedEdge> vec_ed;
  ReducedEdge edge(std::vector<Index>{0, 1, 2});
  vec_ed.push_back(edge);

  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  gn.add(string("name"), string("a"));
  gn1.add(string("name"), string("c"));
  gn2.add(string("name"), string("b"));

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;

  ReducedGraph graph(vec_ed, m_gn);

  Index starting_vertex = 0;
  Branch branch(edge, starting_vertex, graph);

  branch.reverseSequence();
  BOOST_TEST(branch.getBranchStringId() == string("name=b;name=c;name=a;"));
}

BOOST_AUTO_TEST_SUITE_END()
