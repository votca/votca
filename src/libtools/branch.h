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

#include <algorithm>
#include <string>
#include <vector>

#include "votca/tools/graph.h"
#include "votca/tools/reducededge.h"
#include "votca/tools/types.h"

namespace votca {
namespace tools {
class Branch {
 private:
  std::vector<Index> vertex_sequence_;
  std::unordered_map<Index, std::string> node_str_ids_;

 public:
  Branch(const ReducedEdge& edge, const Index starting_vertex,
         const Graph& graph)
      : Branch(edge.getChain(), starting_vertex, graph){};
  // Create a branch, the branches are ordered such that the starting
  // vertex is placed first in the sequence
  Branch(const std::vector<Index>& branch_vertices, const Index starting_vertex,
         const Graph& graph) {

    assert((branch_vertices.front() == starting_vertex ||
            branch_vertices.back() == starting_vertex) &&
           "Cannot create branch with the provided sequence, the provided "
           "starging vertex is not one of the end points");

    vertex_sequence_ = branch_vertices;
    if (branch_vertices.back() == starting_vertex) {
      reverse(vertex_sequence_.begin(), vertex_sequence_.end());
    }
    // Create the string id
    for (Index& vert : vertex_sequence_) {
      GraphNode gn = graph.getNode(vert);
      node_str_ids_[vert] = gn.getStringId();
    }
  }

  void reverseSequence() {
    std::reverse(vertex_sequence_.begin(), vertex_sequence_.end());
  }

  Index getSource() const { return vertex_sequence_.front(); }
  Index getTerminal() const { return vertex_sequence_.back(); }
  std::string getBranchStringId() const noexcept {
    std::string branch_str_id_ = "";
    for (const Index& vert : vertex_sequence_) {
      branch_str_id_ += node_str_ids_.at(vert);
    }
    return branch_str_id_;
  }

  std::vector<Index> getSequence() const noexcept { return vertex_sequence_; }
};

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_BRANCH_H
