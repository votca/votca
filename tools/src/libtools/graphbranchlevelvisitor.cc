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

#include "votca/tools/graphbranchlevelvisitor.h"
#include "votca/tools/graph_bf_visitor.h"
#include "votca/tools/graphvisitor.h"

#include "votca/tools/edge.h"
#include "votca/tools/graph.h"
#include "votca/tools/graphnode.h"

using namespace std;

namespace votca {
namespace tools {

static const string level = "Level";

void GraphBranchLevelVisitor::explore(pair<Index, GraphNode>& p_gn, Graph* g,
                                      const Edge& ed) {
  throw std::runtime_error(
      "GraphBranchLevelVisitor must be able to alter the edge when exploring "
      "it. You have passed in a constant edge");
}

// Add the distance to the node that has not yet been explored
void GraphBranchLevelVisitor::explore(pair<Index, GraphNode>& p_gn, Graph* g,
                                      Edge& ed) {
  // Determine if the node has already been explored
  Index vertex = p_gn.first;
  if (vertex == startingVertex_) {
    ed.add(level, int(0));
  } else {
    // Node has not been explored
    if (explored_.count(vertex) == 0) {
      Index prev_vertex = ed.getOtherEndPoint(vertex);
      // Get the edges associated with this vertex
      std::vector<Edge> prev_edges = g->getNeighEdges(prev_vertex);
      // Cycle them and find the one with the smallest level if that edge
      // even has a level attribute
      Index min_level = -1;
      bool min_level_init = false;
      for (Edge& ed_prev : prev_edges) {
        if (ed_prev != ed) {
          if (ed_prev.exists(level)) {
            if (min_level_init == true) {
              if (ed_prev.get<int>(level) < min_level) {
                min_level = ed_prev.get<int>(level);
              }
            } else {
              min_level = ed.get<int>(level);
              min_level_init = true;
            }
          }
        }
      }  // for

      ed.add(level, min_level++);
    }
  }
  // Ensure the graph node is set to explored
  GraphVisitor::explore(p_gn, g);
}
}  // namespace tools
}  // namespace votca
