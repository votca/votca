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

#include <deque>   // for deque
#include <memory>  // for allocator_traits<>::value_...
#include <queue>   // for queue
#include <set>     // for set
#include <vector>  // for vector

#include "../../include/votca/tools/edge.h"
#include "../../include/votca/tools/graph.h"
#include "../../include/votca/tools/graph_bf_visitor.h"

using namespace std;

namespace votca {
namespace tools {

bool Graph_BF_Visitor::queEmpty() const { return edge_que_.empty(); }

Edge Graph_BF_Visitor::getEdge_(const Graph& graph) {
  Edge oldest_edge = edge_que_.at(0).front();
  edge_que_.at(0).pop();
  if (edge_que_.at(0).size() == 0) {
    edge_que_.pop_front();
  }
  return oldest_edge;
}

// Add edges to be explored
void Graph_BF_Visitor::addEdges_(const Graph& graph, int vertex) {

  vector<Edge> newest_edges = graph.getNeighEdges(vertex);

  // If first edges to be added
  if (edge_que_.empty()) {
    queue<Edge> first_edge_queue;
    for (const Edge edge : newest_edges) {
      int neigh_vert = edge.getOtherEndPoint(vertex);
      if (explored_.count(neigh_vert) == 0) {
        first_edge_queue.push(edge);
      }
    }
    if (!first_edge_queue.empty()) {
      edge_que_.push_back(first_edge_queue);
    }
  } else {

    if (edge_que_.size() == 1) {
      queue<Edge> new_edge_queue;
      for (const Edge edge : newest_edges) {
        int neigh_vert = edge.getOtherEndPoint(vertex);
        if (explored_.count(neigh_vert) == 0) {
          new_edge_queue.push(edge);
        }
      }
      if (!new_edge_queue.empty()) {
        edge_que_.push_back(new_edge_queue);
      }
    } else {
      for (const Edge edge : newest_edges) {
        int neigh_vert = edge.getOtherEndPoint(vertex);
        if (explored_.count(neigh_vert) == 0) {
          edge_que_.at(1).push(edge);
        }
      }
    }
  }
}
}  // namespace tools
}  // namespace votca
