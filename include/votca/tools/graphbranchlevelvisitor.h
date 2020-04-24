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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_TOOLS_GRAPHBRANCHLEVELVISITOR_H
#define VOTCA_TOOLS_GRAPHBRANCHLEVELVISITOR_H
#pragma once

#include "graph_bf_visitor.h"

namespace votca {
namespace tools {

class GraphBranchLevelVisitor : public Graph_BF_Visitor {

 public:
  GraphBranchLevelVisitor() = default;

  /// Note the only manipulation to the BF visitor is the need to add a
  /// level attribute attribute to each of the edges.
  void explore(std::pair<Index, GraphNode>& p_gn, Graph* g, Edge& ed) override;

  void explore(std::pair<Index, GraphNode>& p_gn, Graph* g,
               const Edge& ed = DUMMY_EDGE) override;
};

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_GRAPHBRANCHLEVELVISITOR_H
