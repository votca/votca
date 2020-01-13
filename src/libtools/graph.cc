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

#include "../../include/votca/tools/graph.h"
#include <algorithm>
#include <cassert>
#include <string>

using namespace std;

namespace votca {
namespace tools {

class GraphNode;

bool nodeForEveryVertex_(vector<Index> vertices,
                         unordered_map<Index, GraphNode> nodes) {
  for (const Index vertex : vertices) {
    if (nodes.count(vertex) == 0) {
      return false;
    }
  }
  return true;
}

Graph::Graph(vector<Edge> edges, unordered_map<Index, GraphNode> nodes) {
  edge_container_ = EdgeContainer(edges);
  vector<Index> vertices = edge_container_.getVertices();
  assert(nodeForEveryVertex_(vertices, nodes) &&
         "A node must exist for every vertex.");
  nodes_ = nodes;
  for (const pair<Index, GraphNode> id_and_node : nodes_) {
    if (edge_container_.vertexExist(id_and_node.first) == false) {
      edge_container_.addVertex(id_and_node.first);
    }
  }
  calcId_();
}

bool Graph::operator!=(const Graph& graph) const {
  return id_.compare(graph.id_);
}

bool Graph::operator==(const Graph& graph) const { return !(*(this) != graph); }

vector<pair<Index, GraphNode>> Graph::getIsolatedNodes(void) const {
  vector<pair<Index, GraphNode>> isolated_nodes;
  vector<Index> vertices_degree_0 = edge_container_.getVerticesDegree(0);

  for (const Index vertex : vertices_degree_0) {
    pair<Index, GraphNode> id_and_node{vertex, nodes_.at(vertex)};
    isolated_nodes.push_back(id_and_node);
  }
  return isolated_nodes;
}

vector<pair<Index, GraphNode>> Graph::getNeighNodes(Index vertex) const {
  vector<Index> neigh_vertices = edge_container_.getNeighVertices(vertex);
  vector<pair<Index, GraphNode>> neigh_ids_and_nodes;
  for (const Index& neigh_vert : neigh_vertices) {
    auto id_and_node =
        pair<Index, GraphNode>(neigh_vert, nodes_.at(neigh_vert));
    neigh_ids_and_nodes.push_back(id_and_node);
  }
  return neigh_ids_and_nodes;
}

void Graph::setNode(Index vertex, GraphNode& graph_node) {
  assert(nodes_.count(vertex) && "Can only set a node that already exists");
  nodes_[vertex] = graph_node;
  calcId_();
}

void Graph::setNode(std::pair<Index, GraphNode>& id_and_node) {
  setNode(id_and_node.first, id_and_node.second);
}

GraphNode Graph::getNode(const Index vertex) const {
  assert(nodes_.count(vertex));
  return nodes_.at(vertex);
}

vector<pair<Index, GraphNode>> Graph::getNodes(void) const {
  vector<pair<Index, GraphNode>> vec_nodes;
  for (const pair<Index, GraphNode>& id_and_node : nodes_) {
    vec_nodes.push_back(id_and_node);
  }
  return vec_nodes;
}

vector<Index> Graph::getJunctions() const {
  vector<Index> junctions;
  Index max_degree = edge_container_.getMaxDegree();
  for (Index degree = 3; degree <= max_degree; ++degree) {
    vector<Index> vertices = edge_container_.getVerticesDegree(degree);
    junctions.insert(junctions.end(), vertices.begin(), vertices.end());
  }
  return junctions;
}

void Graph::clearNodes() { nodes_.clear(); }

void Graph::copyNodes(Graph& graph) {
  assert(nodes_.size() == 0);
  for (const pair<Index, GraphNode>& id_and_node : graph.nodes_) {
    this->nodes_[id_and_node.first] = id_and_node.second;
  }
}

void Graph::calcId_() {
  auto nodes = getNodes();
  sort(nodes.begin(), nodes.end(), cmpVertNodePair);
  string struct_Id_temp = "";
  for (const pair<Index, GraphNode>& id_and_node : nodes) {
    struct_Id_temp.append(id_and_node.second.getStringId());
  }
  id_ = struct_Id_temp;
  return;
}

Index Graph::getDegree(Index vertex) const {
  if (edge_container_.vertexExist(vertex)) {
    return edge_container_.getDegree(vertex);
  }
  if (nodes_.count(vertex)) {
    return 0;
  }
  throw invalid_argument(
      "vertex does not exist within the graph the degree is "
      "not defined.");
}

bool Graph::vertexExist(Index vertex) const {
  if (edge_container_.vertexExist(vertex)) {
    return true;
  }
  if (nodes_.count(vertex)) {
    return true;
  }
  return false;
}

vector<Index> Graph::getVerticesDegree(Index degree) const {
  return edge_container_.getVerticesDegree(degree);
}

vector<Index> Graph::getVertices() const {
  return edge_container_.getVertices();
}

ostream& operator<<(ostream& os, const Graph graph) {
  os << "Graph" << endl;
  for (const pair<Index, GraphNode>& id_and_node : graph.nodes_) {
    os << "Node " << id_and_node.first << endl;
    os << id_and_node.second << endl;
  }
  return os;
}

bool cmpVertNodePair(const pair<Index, GraphNode>& id_and_node1,
                     const pair<Index, GraphNode>& id_and_node2) {
  string str1_Id = id_and_node1.second.getStringId();
  return str1_Id.compare(id_and_node2.second.getStringId()) < 0;
}
}  // namespace tools
}  // namespace votca
