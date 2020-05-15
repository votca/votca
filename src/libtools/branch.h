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
 * Unless required by applicable law or agreed to in writingraph, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_TOOLS_BRANCH_H
#define VOTCA_TOOLS_BRANCH_H
#pragma once

#include <string>
#include <vector>

#include "votca/tools/graph.h"
#include "votca/tools/reducededge.h"
#include "votca/tools/types.h"

namespace votca {
namespace tools {

class Level {
 public:
  int num;
  Level() : num(-1){};
  Level(int lev_num) : num(lev_num){};

  bool operator==(const Level& lev) const { return num == lev.num; }
  bool operator!=(const Level& lev) const { return !(*this == lev); }
  bool operator<(const Level& lev) const { return num < lev.num; }
  bool operator>(const Level& lev) const {
    return !(*this < lev) && (*this != lev);
  }
  bool operator<=(const Level& lev) const { return !(*this > lev); }
  bool operator>=(const Level& lev) const { return !(*this < lev); }
};

}  // namespace tools
}  // namespace votca
// custom specialization of std::hash can be injected in namespace std
namespace std {
template <>
struct hash<votca::tools::Level> {
  std::size_t operator()(votca::tools::Level const& level) const noexcept {
    return std::hash<int>{}(level.num);
  }
};
}  // namespace std

namespace votca {
namespace tools {
class Branch {
 private:
  std::vector<Index> vertex_sequence_;
  std::unordered_map<Index, ContentLabel> node_labels_;
  Level level_;

  void init_(const std::vector<Index>& branch_vertices,
             const Index starting_vertex, const Graph& graph) {
    assert((branch_vertices.front() == starting_vertex ||
            branch_vertices.back() == starting_vertex) &&
           "Cannot create branch with the provided sequence, the provided "
           "starging vertex is not one of the end points");

    vertex_sequence_ = branch_vertices;
    if (branch_vertices.back() == starting_vertex) {
      reverse(vertex_sequence_.begin(), vertex_sequence_.end());
    }
    // Add the content label of the node
    for (Index& vert : vertex_sequence_) {
      GraphNode gn = graph.getNode(vert);
      node_labels_[vert] = gn.getContentLabel();
    }
  }

 public:
  Branch(const ReducedEdge& edge, const Index starting_vertex,
         const Graph& graph) {

    if (not edge.exists("Level")) {
      throw std::runtime_error(
          "Cannot create branch from reduced edge it has not been assigned a "
          "level");
    }
    level_ = Level(edge.get<int>("Level"));
    init_(edge.getChain(), starting_vertex, graph);
  };
  // Create a branch, the branches are ordered such that the starting
  // vertex is placed first in the sequence
  Branch(const std::vector<Index>& branch_vertices, int level,
         const Index starting_vertex, const Graph& graph) {

    level_ = Level(level);
    init_(branch_vertices, starting_vertex, graph);
  }

  Level getLevel() const noexcept { return level_; }

  void reverseSequence() {
    std::reverse(vertex_sequence_.begin(), vertex_sequence_.end());
  }

  Index getSource() const { return vertex_sequence_.front(); }
  Index getTerminal() const { return vertex_sequence_.back(); }

  ContentLabel getContentLabel() const noexcept {
    ContentLabel label;
    std::cout << "getting content label from branch" << std::endl;
    for (const Index& vert : vertex_sequence_) {
      label.append(node_labels_.at(vert));
    }
    std::cout << "Making label branch label" << std::endl;
    label.makeBranch();
    return label;
  }

  std::vector<Index> getSequence() const noexcept { return vertex_sequence_; }
};

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_BRANCH_H
