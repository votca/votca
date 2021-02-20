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

#define BOOST_TEST_MODULE branchsequenceleaf_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

#include "../libtools/branch.h"
#include "../libtools/branchsequenceleaf.h"
#include "votca/tools/graphnode.h"
#include "votca/tools/reducededge.h"
#include "votca/tools/reducedgraph.h"

using namespace std;
using namespace votca;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(branchsequenceleaf_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  vector<ReducedEdge> vec_ed;
  ReducedEdge edge(std::vector<Index>{0, 1, 2});
  edge.add("Level", 0);
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
  std::cout << "Getting branch content label string " << std::endl;
  std::cout << branch.getContentLabel().get() << std::endl;
  Index branch_id = 0;
  BranchSequenceLeaf sequenceleaf(branch_id, branch);
}

BOOST_AUTO_TEST_CASE(canOrder_isDangling_makeBranchEnd) {
  vector<ReducedEdge> vec_ed;
  ReducedEdge edge(std::vector<Index>{0, 1, 2});
  edge.add("Level", 0);
  vec_ed.push_back(edge);

  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;

  ReducedGraph graph(vec_ed, m_gn);

  Index starting_vertex = 0;
  Branch branch(edge, starting_vertex, graph);

  Index branch_id = 0;
  BranchSequenceLeaf sequenceleaf(branch_id, branch);

  BOOST_TEST(sequenceleaf.canOrder() == false);
  BOOST_TEST(sequenceleaf.isDangling() == true);

  sequenceleaf.makeBranchEnd();

  BOOST_TEST(sequenceleaf.canOrder() == true);
  BOOST_TEST(sequenceleaf.isDangling() == false);
}

BOOST_AUTO_TEST_CASE(canAddLeaf_sequenceIsIncomplete_addLeaf) {
  vector<ReducedEdge> vec_ed;
  /// First branch
  ReducedEdge edge(std::vector<Index>{0, 1, 2});
  edge.add("Level", 0);
  /// Second branch
  ReducedEdge edge1(std::vector<Index>{2, 3});
  edge1.add("Level", 0);
  /// Third branch
  ReducedEdge edge2(std::vector<Index>{2, 4, 5});
  edge2.add("Level", 0);

  vec_ed.push_back(edge);
  vec_ed.push_back(edge1);
  vec_ed.push_back(edge2);

  /// First branch
  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  /// Second branch
  GraphNode gn3;

  /// Fourth branch
  GraphNode gn4;
  GraphNode gn5;

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;
  m_gn[3] = gn3;
  m_gn[4] = gn4;
  m_gn[5] = gn5;

  ReducedGraph graph(vec_ed, m_gn);

  Index starting_vertex = 0;
  Branch branch(edge, starting_vertex, graph);
  starting_vertex = 2;
  Branch branch1(edge1, starting_vertex, graph);
  starting_vertex = 5;
  Branch branch2(edge2, starting_vertex, graph);

  Index branch_id = 0;
  BranchSequenceLeaf sequenceleaf(branch_id, branch);

  Index branch_id1 = 1;
  BOOST_TEST(sequenceleaf.sequenceIsIncomplete() == true);
  auto leaf1 = std::make_shared<BranchSequenceLeaf>(branch_id1, branch1);
  BOOST_TEST(sequenceleaf.canAddLeaf(leaf1) == true);

  sequenceleaf.addLeaf(leaf1);

  Index branch_id2 = 2;
  BOOST_TEST(sequenceleaf.sequenceIsIncomplete() == true);
  auto leaf2 = std::make_shared<BranchSequenceLeaf>(branch_id2, branch2);
  BOOST_TEST(sequenceleaf.canAddLeaf(leaf2) == true);

  sequenceleaf.addLeaf(leaf2);

  BOOST_TEST(sequenceleaf.sequenceIsIncomplete() == false);
  BOOST_TEST(sequenceleaf.treeIsIncomplete() == true);
}

BOOST_AUTO_TEST_CASE(getTreeContentLabel) {
  vector<ReducedEdge> vec_ed;
  /// First branch
  ReducedEdge edge(std::vector<Index>{0, 1, 2});
  edge.add("Level", 1);
  /// Second branch
  ReducedEdge edge1(std::vector<Index>{2, 3});
  edge1.add("Level", 1);
  /// Third branch
  ReducedEdge edge2(std::vector<Index>{2, 4, 5});
  edge2.add("Level", 1);

  vec_ed.push_back(edge);
  vec_ed.push_back(edge1);
  vec_ed.push_back(edge2);

  /// First branch
  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  /// Second branch
  GraphNode gn3;

  /// Fourth branch
  GraphNode gn4;
  GraphNode gn5;

  gn.add(string("name"), string("a"));
  gn1.add(string("name"), string("a"));
  gn2.add(string("name"), string("a"));
  gn3.add(string("name"), string("c"));
  gn4.add(string("name"), string("b"));
  gn5.add(string("name"), string("b"));

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;
  m_gn[3] = gn3;
  m_gn[4] = gn4;
  m_gn[5] = gn5;

  {  // First test
    ReducedGraph graph(vec_ed, m_gn);

    Index starting_vertex = 0;
    Branch branch(edge, starting_vertex, graph);
    starting_vertex = 2;
    Branch branch1(edge1, starting_vertex, graph);
    starting_vertex = 5;
    Branch branch2(edge2, starting_vertex, graph);

    Index branch_id = 0;
    BranchSequenceLeaf sequenceleaf(branch_id, branch);

    Index branch_id1 = 1;
    auto leaf1 = std::make_shared<BranchSequenceLeaf>(branch_id1, branch1);
    sequenceleaf.addLeaf(leaf1);
    Index branch_id2 = 2;
    auto leaf2 = std::make_shared<BranchSequenceLeaf>(branch_id2, branch2);
    sequenceleaf.addLeaf(leaf2);

    sequenceleaf.sortBranchSequence();

    std::string test_str =
        "{name=a;name=a;name=a}{name=a;name=c}{name=a;name=b;name=b}";
    std::string cum_id = sequenceleaf.getTreeContentLabel().get();
    std::cout << test_str << std::endl;
    std::cout << cum_id << std::endl;
    BOOST_TEST(cum_id.compare(test_str) == 0);
  }
  /// Another graph node to the second branch so that they are equal length
  /// Second branch
  GraphNode gn6;

  gn6.add(string("name"), string("c"));
  m_gn[6] = gn6;

  {  // Second test
    ReducedGraph graph(vec_ed, m_gn);

    // Replace edge1 to include vertex 6
    edge1 = ReducedEdge(std::vector<Index>{2, 3, 6});
    edge1.add("Level", 1);
    Index starting_vertex = 0;
    Branch branch(edge, starting_vertex, graph);
    starting_vertex = 2;
    Branch branch1(edge1, starting_vertex, graph);
    starting_vertex = 5;
    Branch branch2(edge2, starting_vertex, graph);

    Index branch_id = 0;
    BranchSequenceLeaf sequenceleaf(branch_id, branch);

    Index branch_id1 = 1;
    auto leaf1 = std::make_shared<BranchSequenceLeaf>(branch_id1, branch1);
    sequenceleaf.addLeaf(leaf1);
    Index branch_id2 = 2;
    auto leaf2 = std::make_shared<BranchSequenceLeaf>(branch_id2, branch2);
    sequenceleaf.addLeaf(leaf2);

    sequenceleaf.sortBranchSequence();
    std::string cum_id = sequenceleaf.getTreeContentLabel().get();
    BOOST_TEST(cum_id == std::string("{name=a;name=a;name=a}{name=a;name=b;"
                                     "name=b}{name=a;name=c;name=c}"));
  }
}

BOOST_AUTO_TEST_CASE(getCanonicalVertexSequence) {
  vector<ReducedEdge> vec_ed;
  /// First branch
  ReducedEdge edge(std::vector<Index>{0, 1, 2});
  edge.add("Level", 1);
  /// Second branch
  ReducedEdge edge1(std::vector<Index>{2, 3});
  edge1.add("Level", 1);
  /// Third branch
  ReducedEdge edge2(std::vector<Index>{2, 4, 5});
  edge2.add("Level", 1);

  vec_ed.push_back(edge);
  vec_ed.push_back(edge1);
  vec_ed.push_back(edge2);

  /// First branch
  GraphNode gn;
  GraphNode gn1;
  GraphNode gn2;

  /// Second branch
  GraphNode gn3;

  /// Fourth branch
  GraphNode gn4;
  GraphNode gn5;

  gn.add(string("name"), string("a"));
  gn1.add(string("name"), string("a"));
  gn2.add(string("name"), string("a"));
  gn3.add(string("name"), string("c"));
  gn4.add(string("name"), string("b"));
  gn5.add(string("name"), string("b"));

  unordered_map<votca::Index, GraphNode> m_gn;
  m_gn[0] = gn;
  m_gn[1] = gn1;
  m_gn[2] = gn2;
  m_gn[3] = gn3;
  m_gn[4] = gn4;
  m_gn[5] = gn5;

  ///
  ///    branch 0     branch 2
  /// 0     1     2     4     5
  /// a  -  a  -  a  -  b  -  b
  ///             |
  ///             c 3
  ///         branch 1
  ///

  {  // First test
    ReducedGraph graph(vec_ed, m_gn);

    Index starting_vertex = 0;
    Branch branch(edge, starting_vertex, graph);
    starting_vertex = 2;
    Branch branch1(edge1, starting_vertex, graph);
    starting_vertex = 5;
    Branch branch2(edge2, starting_vertex, graph);

    Index branch_id = 0;
    BranchSequenceLeaf sequenceleaf(branch_id, branch);

    Index branch_id1 = 1;
    auto leaf1 = std::make_shared<BranchSequenceLeaf>(branch_id1, branch1);
    sequenceleaf.addLeaf(leaf1);
    Index branch_id2 = 2;
    auto leaf2 = std::make_shared<BranchSequenceLeaf>(branch_id2, branch2);
    sequenceleaf.addLeaf(leaf2);

    sequenceleaf.sortBranchSequence();
    std::vector<Index> seq = sequenceleaf.getCanonicalVertexSequence();

    std::vector<Index> correct_sequence = {0, 1, 2, 3, 4, 5};
    BOOST_TEST(correct_sequence.size() == seq.size());
    for (size_t i = 0; i < correct_sequence.size(); ++i) {
      BOOST_TEST(correct_sequence.at(i) == seq.at(i));
    }
  }
  /// Another graph node to the second branch so that they are equal length
  /// Second branch
  GraphNode gn6;

  gn6.add(string("name"), string("c"));
  m_gn[6] = gn6;

  ///
  ///    branch 0     branch 2
  /// 0     1     2     4     5
  /// a  -  a  -  a  -  b  -  b
  ///             |
  ///             c 3
  ///             |
  ///             c 6
  ///         branch 1
  ///

  {  // Second test
    ReducedGraph graph(vec_ed, m_gn);

    // Replace edge1 to include vertex 6
    edge1 = ReducedEdge(std::vector<Index>{2, 3, 6});
    edge1.add("Level", 1);
    Index starting_vertex = 0;
    Branch branch(edge, starting_vertex, graph);
    starting_vertex = 2;
    Branch branch1(edge1, starting_vertex, graph);
    /// It should not matter that 5 is set as the starting vertex, as long
    /// as 5 is one of the ends
    starting_vertex = 5;
    Branch branch2(edge2, starting_vertex, graph);

    Index branch_id = 0;
    BranchSequenceLeaf sequenceleaf(branch_id, branch);

    Index branch_id1 = 1;
    auto leaf1 = std::make_shared<BranchSequenceLeaf>(branch_id1, branch1);
    sequenceleaf.addLeaf(leaf1);
    Index branch_id2 = 2;
    auto leaf2 = std::make_shared<BranchSequenceLeaf>(branch_id2, branch2);
    sequenceleaf.addLeaf(leaf2);

    sequenceleaf.sortBranchSequence();
    std::vector<Index> seq = sequenceleaf.getCanonicalVertexSequence();

    std::vector<Index> correct_sequence = {0, 1, 2, 4, 5, 3, 6};
    BOOST_TEST(correct_sequence.size() == seq.size());
    for (size_t i = 0; i < correct_sequence.size(); ++i) {
      BOOST_TEST(correct_sequence.at(i) == seq.at(i));
    }
  }

  ///
  ///    branch 0       -  b 4 -
  /// 0     1     2  /           \
  /// a  -  a  -  a     branch 2
  ///             |  \           /
  ///             c 3   -  b 5 -
  ///             |
  ///             c 6
  ///         branch 1
  ///

  {  // third test
    ReducedGraph graph(vec_ed, m_gn);

    // Replace edge1 to include vertex 6
    edge1 = ReducedEdge(std::vector<Index>{2, 3, 6});
    edge1.add("Level", 1);
    edge2 = ReducedEdge(std::vector<Index>{2, 4, 5, 2});
    edge2.add("Level", 1);

    Index starting_vertex = 0;
    Branch branch(edge, starting_vertex, graph);
    starting_vertex = 2;
    Branch branch1(edge1, starting_vertex, graph);
    starting_vertex = 2;
    Branch branch2(edge2, starting_vertex, graph);

    Index branch_id = 0;
    BranchSequenceLeaf sequenceleaf(branch_id, branch);

    Index branch_id1 = 1;
    auto leaf1 = std::make_shared<BranchSequenceLeaf>(branch_id1, branch1);
    sequenceleaf.addLeaf(leaf1);
    Index branch_id2 = 2;
    auto leaf2 = std::make_shared<BranchSequenceLeaf>(branch_id2, branch2);
    sequenceleaf.addLeaf(leaf2);

    sequenceleaf.sortBranchSequence();
    std::vector<Index> seq = sequenceleaf.getCanonicalVertexSequence();

    // Here the branch 2 should come come last because it is now longer
    // than the branch 1 because the starting vertex is included twice,
    // also the order of vertices 4 and 5 is irrelevant

    BOOST_TEST(seq.at(0) == 0);
    BOOST_TEST(seq.at(1) == 1);
    BOOST_TEST(seq.at(2) == 2);
    BOOST_TEST(seq.at(3) == 3);
    BOOST_TEST(seq.at(4) == 6);
    BOOST_TEST((seq.at(5) == 5 || seq.at(5) == 4));
    BOOST_TEST((seq.at(6) == 4 || seq.at(6) == 5));
  }

  /// Another graph node to the first branch so that they are equal length
  /// Second branch
  GraphNode gn7;

  gn7.add(string("name"), string("a"));
  m_gn[7] = gn7;

  ///
  ///    branch 0       -  b 4 -
  /// 0     1     2  /           \
    /// a  -  a  -  a     branch 2
  ///             |  \           /
  ///             c 3   -  b 5 -
  ///             |
  ///             c 6
  ///             |
  ///             a 7
  ///         branch 1
  ///

  {  // third test
    ReducedGraph graph(vec_ed, m_gn);

    // Replace edge1 to include vertex 6
    edge1 = ReducedEdge(std::vector<Index>{2, 3, 6, 7});
    edge1.add("Level", 1);
    edge2 = ReducedEdge(std::vector<Index>{2, 4, 5, 2});
    edge2.add("Level", 1);

    Index starting_vertex = 0;
    Branch branch(edge, starting_vertex, graph);
    starting_vertex = 2;
    Branch branch1(edge1, starting_vertex, graph);
    starting_vertex = 2;
    Branch branch2(edge2, starting_vertex, graph);

    Index branch_id = 0;
    BranchSequenceLeaf sequenceleaf(branch_id, branch);

    Index branch_id1 = 1;
    auto leaf1 = std::make_shared<BranchSequenceLeaf>(branch_id1, branch1);
    sequenceleaf.addLeaf(leaf1);
    Index branch_id2 = 2;
    auto leaf2 = std::make_shared<BranchSequenceLeaf>(branch_id2, branch2);
    sequenceleaf.addLeaf(leaf2);

    sequenceleaf.sortBranchSequence();
    std::vector<Index> seq = sequenceleaf.getCanonicalVertexSequence();

    for (Index val : seq) {
      std::cout << "val is " << val << std::endl;
    }
    // Here the branch 2 should come come last because it is now longer
    // than the branch 1 because the starting vertex is included twice,
    // also the order of vertices 4 and 5 is irrelevant

    BOOST_TEST(seq.at(0) == 0);
    BOOST_TEST(seq.at(1) == 1);
    BOOST_TEST(seq.at(2) == 2);
    BOOST_TEST((seq.at(3) == 5 || seq.at(3) == 4));
    BOOST_TEST((seq.at(4) == 4 || seq.at(4) == 5));
    BOOST_TEST(seq.at(5) == 3);
    BOOST_TEST(seq.at(6) == 6);
    BOOST_TEST(seq.at(7) == 7);
  }

  ///
  ///              a 0
  ///               |
  ///      - - - - a 1 - - - -
  ///     |         |         |
  ///    a 2 - - - a 3 - - - a 4
  ///     |         |         |
  ///      - - - - a 5 - - - -
  ///
  /// branch 0  a0, a1
  /// branch 1  a1, a2
  /// branch 2  a1, a3
  /// branch 3  a1, a4
  /// branch 4  a2, a3
  /// branch 5  a4, a3
  /// branch 6  a2, a5
  /// branch 7  a3, a5
  /// branch 8  a4, a5
  ///
  /// The purpose of this branch is to distinguish between
  /// Branches 1, 2 and 3, branch 2 is different from branches 1 and 3 and
  /// this should show up in the vertex sequence
  ///

  vector<ReducedEdge> vec_ed2;
  /// First branch
  edge = ReducedEdge(std::vector<Index>{0, 1});
  edge.add("Level", 1);
  /// Second branch
  edge1 = ReducedEdge(std::vector<Index>{1, 2});
  edge1.add("Level", 1);
  /// Third branch
  edge2 = ReducedEdge(std::vector<Index>{1, 3});
  edge2.add("Level", 1);
  /// Fourth branch
  ReducedEdge edge3 = ReducedEdge(std::vector<Index>{1, 4});
  edge3.add("Level", 1);
  /// Fith branch
  ReducedEdge edge4 = ReducedEdge(std::vector<Index>{2, 3});
  edge4.add("Level", 1);
  /// Sixth branch
  ReducedEdge edge5 = ReducedEdge(std::vector<Index>{4, 3});
  edge5.add("Level", 1);
  /// Seventh branch
  ReducedEdge edge6 = ReducedEdge(std::vector<Index>{2, 5});
  edge6.add("Level", 1);
  /// eigth branch
  ReducedEdge edge7 = ReducedEdge(std::vector<Index>{3, 5});
  edge7.add("Level", 1);
  /// ninth branch
  ReducedEdge edge8 = ReducedEdge(std::vector<Index>{4, 5});
  edge8.add("Level", 1);

  vec_ed.push_back(edge);
  vec_ed.push_back(edge1);
  vec_ed.push_back(edge2);
  vec_ed.push_back(edge3);
  vec_ed.push_back(edge4);
  vec_ed.push_back(edge5);
  vec_ed.push_back(edge6);
  vec_ed.push_back(edge7);
  vec_ed.push_back(edge8);

  gn3.reset(std::string("name"), std::string("a"));
  gn4.reset(std::string("name"), std::string("a"));
  gn5.reset(std::string("name"), std::string("a"));

  unordered_map<votca::Index, GraphNode> m_gn2;
  m_gn2[0] = gn;
  m_gn2[1] = gn1;
  m_gn2[2] = gn2;
  m_gn2[3] = gn3;
  m_gn2[4] = gn4;
  m_gn2[5] = gn5;

  {
    ReducedGraph graph(vec_ed, m_gn);

    // Start by creating branch sequences from branches 6, 7 and 8 and making
    // them complete, end nodes
    ///
    ///              a 0
    ///               |
    ///      - - - - a 1 - - - -
    ///     |         |         |
    ///    a 2 - - - a 3 - - - a 4
    ///     |         |         |
    ///      - - - - a 5 - - - -
    ///
    /// branch 0  a0, a1  branch dist 0
    /// branch 1  a1, a2  branch dist 1
    /// branch 2  a1, a3  branch dist 1
    /// branch 3  a1, a4  branch dist 1
    /// branch 4  a2, a3  branch dist 2
    /// branch 5  a4, a3  branch dist 2
    /// branch 6  a2, a5  branch dist 2
    /// branch 7  a3, a5  branch dist 2
    /// branch 8  a4, a5  branch dist 2

    Index starting_vertex = 2;
    Branch branch6(edge6, starting_vertex, graph);
    starting_vertex = 3;
    Branch branch7(edge7, starting_vertex, graph);
    starting_vertex = 4;
    Branch branch8(edge8, starting_vertex, graph);

    Index branch_id = 6;
    BranchSequenceLeaf sequenceleaf6(branch_id, branch6);
    branch_id = 7;
    BranchSequenceLeaf sequenceleaf7(branch_id, branch7);
    branch_id = 8;
    BranchSequenceLeaf sequenceleaf8(branch_id, branch8);

    // in what order do we add the sequence nodes and how should they be
    // linked
  }
}
BOOST_AUTO_TEST_SUITE_END()
