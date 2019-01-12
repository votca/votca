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

#include <votca/tools/graph_bf_visitor.h>
#include <votca/tools/graphdistvisitor.h>
#include <votca/tools/graphvisitor.h>

#include <votca/tools/edge.h>
#include <votca/tools/graph.h>
#include <votca/tools/graphnode.h>

using namespace std;

namespace votca {
namespace tools {

// Add the distance to the node that has not yet been explored
void GraphDistVisitor::exploreNode(pair<int, GraphNode> &p_gn, Graph& g,
                                    Edge ed) {
  // Determine if the node has already been explored
  int vertex = p_gn.first;
  if (vertex == startingVertex_) {
    p_gn.second.int_vals_["Dist"] = 0;
    p_gn.second.initStringId_();
    // Update the graph with new graph node
    g.setNode(p_gn);
  } else {
    // Node has not been explored
    if (explored_.count(vertex) == 0) {
      int prev_vertex = ed.getOtherEndPoint(vertex);
      GraphNode gn_prev = g.getNode(prev_vertex);
      p_gn.second.int_vals_["Dist"] = gn_prev.int_vals_["Dist"] + 1;
      p_gn.second.initStringId_();
      g.setNode(p_gn);
    }
  }
  // Ensure the graph node is set to explored
  GraphVisitor::exploreNode(p_gn, g);
}
}
}
