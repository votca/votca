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

#include <exception>
#include <iostream>
#include <vector>
#include <votca/tools/edge.h>
#include <votca/tools/graph.h>
#include <votca/tools/graphvisitor.h>

using namespace std;

namespace votca {
namespace tools {

class GraphNode;

bool GraphVisitor::queEmpty() const { return true; }

void GraphVisitor::addEdges_(const Graph& graph, int vertex) {
  throw runtime_error("addEdges_ method must be defined by your visitor");
}

void GraphVisitor::exploreNode(pair<int, GraphNode>& vertex_and_node,
                               Graph& graph, Edge edge) {
  explored_.insert(vertex_and_node.first);
}

vector<int> GraphVisitor::getUnexploredVertex(const Edge edge) const {
  vector<int> unexp_vert;
  if (explored_.count(edge.getEndPoint1()) == 0) {
    unexp_vert.push_back(edge.getEndPoint1());
  }
  if (explored_.count(edge.getEndPoint2()) == 0) {
    unexp_vert.push_back(edge.getEndPoint2());
  }
  return unexp_vert;
}

bool GraphVisitor::vertexExplored(const int vertex) const {
  return explored_.count(vertex) == 1;
}

void GraphVisitor::initialize(Graph& graph) {
  vector<Edge> neigh_eds = graph.getNeighEdges(startingVertex_);
  GraphNode graph_node = graph.getNode(startingVertex_);
  pair<int, GraphNode> vertex_and_graph_node(startingVertex_, graph_node);
  exploreNode(vertex_and_graph_node, graph);
  addEdges_(graph, startingVertex_);
}

void GraphVisitor::exec(Graph& graph, Edge edge) {
  vector<int> unexp_vert = getUnexploredVertex(edge);
  // If no vertices are return than just ignore it means the same
  // vertex was explored from a different direction
  if (!unexp_vert.size()) return;
  // If two values are returned this is a problem
  if (unexp_vert.size() > 1) {
    throw runtime_error(
        "More than one unexplored vertex in an edge,"
        " did you set the starting node");
  }

  pair<int, GraphNode> vertex_and_node(unexp_vert.at(0),
                                       graph.getNode(unexp_vert.at(0)));

  exploreNode(vertex_and_node, graph, edge);
}

Edge GraphVisitor::getEdge_(const Graph& graph) {
  throw runtime_error("The getEdge_ function must be set");
}

Edge GraphVisitor::nextEdge(Graph graph) {

  // Get the edge and at the same time remove it from whatever queue it is in
  Edge edge = getEdge_(graph);
  vector<int> unexplored_vertices = getUnexploredVertex(edge);
  // Do not add neighboring edges if they belong to a vertex that has already
  // been explored because they will have already been added
  if (unexplored_vertices.size()) {
    addEdges_(graph, unexplored_vertices.at(0));
  }
  return edge;
}

set<int> GraphVisitor::getExploredVertices() const { return explored_; }

}  // namespace tools
}  // namespace votca
