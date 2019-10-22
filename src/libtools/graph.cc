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

#include <algorithm>
#include <cassert>
#include <string>
#include <votca/tools/graph.h>

using namespace std;

namespace votca {
namespace tools {

class GraphNode;

bool nodeForEveryVertex_(vector<int> vertices,
                         unordered_map<int, GraphNode> nodes) {
  for (const int vertex : vertices) {
    if (nodes.count(vertex) == 0) {
      return false;
    }
  }
  return true;
}

Graph::Graph(vector<Edge> edges, unordered_map<long int, GraphNode> nodes) {
  edge_container_ = EdgeContainer(edges);
  vector<int> vertices = edge_container_.getVertices();
  assert(nodeForEveryVertex_(vertices, nodes) &&
         "A node must exist for every vertex.");
  nodes_ = nodes;
  for (const pair<int, GraphNode> id_and_node : nodes_) {
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

vector<pair<int, GraphNode>> Graph::getIsolatedNodes(void) const {
  vector<pair<int, GraphNode>> isolated_nodes;
  vector<int> vertices_degree_0 = edge_container_.getVerticesDegree(0);

  for (const int vertex : vertices_degree_0) {
    pair<int, GraphNode> id_and_node{vertex, nodes_.at(vertex)};
    isolated_nodes.push_back(id_and_node);
  }
  return isolated_nodes;
}

vector<pair<int, GraphNode>> Graph::getNeighNodes(int vertex) const {
  vector<int> neigh_vertices = edge_container_.getNeighVertices(vertex);
  vector<pair<int, GraphNode>> neigh_ids_and_nodes;
  for (const int& neigh_vert : neigh_vertices) {
    auto id_and_node = pair<int, GraphNode>(neigh_vert, nodes_.at(neigh_vert));
    neigh_ids_and_nodes.push_back(id_and_node);
  }
  return neigh_ids_and_nodes;
}

void Graph::setNode(int vertex, GraphNode& graph_node) {
  assert(nodes_.count(vertex) && "Can only set a node that already exists");
  nodes_[vertex] = graph_node;
  calcId_();
}

void Graph::setNode(std::pair<int, GraphNode>& id_and_node) {
  setNode(id_and_node.first, id_and_node.second);
}

GraphNode Graph::getNode(const int vertex) const {
  assert(nodes_.count(vertex));
  return nodes_.at(vertex);
}

vector<pair<int, GraphNode>> Graph::getNodes(void) const {
  vector<pair<int, GraphNode>> vec_nodes;
  for (const pair<int, GraphNode>& id_and_node : nodes_) {
    vec_nodes.push_back(id_and_node);
  }
  return vec_nodes;
}

vector<int> Graph::getJunctions() const {
  vector<int> junctions;
  int max_degree = edge_container_.getMaxDegree();
  for (int degree = 3; degree <= max_degree; ++degree) {
    vector<int> vertices = edge_container_.getVerticesDegree(degree);
    junctions.insert(junctions.end(), vertices.begin(), vertices.end());
  }
  return junctions;
}

void Graph::clearNodes() { nodes_.clear(); }

void Graph::copyNodes(Graph& graph) {
  assert(nodes_.size() == 0);
  for (const pair<int, GraphNode>& id_and_node : graph.nodes_) {
    this->nodes_[id_and_node.first] = id_and_node.second;
  }
}

void Graph::calcId_() {
  auto nodes = getNodes();
  sort(nodes.begin(), nodes.end(), cmpVertNodePair);
  string struct_Id_temp = "";
  for (const pair<int, GraphNode>& id_and_node : nodes) {
    struct_Id_temp.append(id_and_node.second.getStringId());
  }
  id_ = struct_Id_temp;
  return;
}

int Graph::getDegree(int vertex) const {
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

bool Graph::vertexExist(int vertex) const {
  if (edge_container_.vertexExist(vertex)) {
    return true;
  }
  if (nodes_.count(vertex)) {
    return true;
  }
  return false;
}

vector<int> Graph::getVerticesDegree(int degree) const {
  return edge_container_.getVerticesDegree(degree);
}

vector<int> Graph::getVertices() const { return edge_container_.getVertices(); }

ostream& operator<<(ostream& os, const Graph graph) {
  os << "Graph" << endl;
  for (const pair<int, GraphNode>& id_and_node : graph.nodes_) {
    os << "Node " << id_and_node.first << endl;
    os << id_and_node.second << endl;
  }
  return os;
}

bool cmpVertNodePair(const pair<int, GraphNode>& id_and_node1,
                     const pair<int, GraphNode>& id_and_node2) {
  string str1_Id = id_and_node1.second.getStringId();
  return str1_Id.compare(id_and_node2.second.getStringId()) < 0;
}
}  // namespace tools
}  // namespace votca
