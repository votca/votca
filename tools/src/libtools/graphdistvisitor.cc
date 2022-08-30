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

// Local VOTCA includes
#include "votca/tools/graphdistvisitor.h"
#include "votca/tools/edge.h"
#include "votca/tools/graph.h"
#include "votca/tools/graph_bf_visitor.h"
#include "votca/tools/graphnode.h"
#include "votca/tools/graphvisitor.h"

using namespace std;

namespace votca {
namespace tools {

static const string dist = "Dist";
// Add the distance to the node that has not yet been explored
void GraphDistVisitor::explore(pair<Index, GraphNode>& p_gn, Graph* g,
                               const Edge& ed) {
  // Determine if the node has already been explored
  Index vertex = p_gn.first;
  if (vertex == startingVertex_) {
    if (p_gn.second.exists(dist)) {
      p_gn.second.reset(dist, int(0));
    } else {
      p_gn.second.add(dist, int(0));
    }
    // Update the graph with new graph node
    g->setNode(p_gn);
  } else {
    // Node has not been explored
    if (explored_.count(vertex) == 0) {
      Index prev_vertex = ed.getOtherEndPoint(vertex);
      GraphNode gn_prev = g->getNode(prev_vertex);
      if (p_gn.second.exists(dist)) {
        p_gn.second.reset(dist, gn_prev.get<int>(dist) + 1);
      } else {
        p_gn.second.add(dist, gn_prev.get<int>(dist) + 1);
      }
      g->setNode(p_gn);
    }
  }
  // Ensure the graph node is set to explored
  GraphVisitor::explore(p_gn, g);
}
}  // namespace tools
}  // namespace votca
