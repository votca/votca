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

#include <algorithm>
#include <cassert>
#include <string>
#include <votca/tools/graph.h>

using namespace std;

namespace votca {
namespace tools {

class GraphNode;

bool Graph::operator!=(const Graph& graph) const {
  return id_.compare(graph.id_);
}

bool Graph::operator==(const Graph& graph) const { return !(*(this) != graph); }

Graph::Graph(const Graph& graph) {
  this->edge_container_ = graph.edge_container_;
  for (const pair<int, GraphNode>& id_and_node : graph.nodes_) {
    this->nodes_[id_and_node.first] = id_and_node.second;
  }
  this->id_ = graph.id_;
}

Graph& Graph::operator=(const Graph& graph) {
  this->edge_container_ = graph.edge_container_;
  for (const pair<int, GraphNode>& id_and_node : graph.nodes_) {
    this->nodes_[id_and_node.first] = id_and_node.second;
  }
  this->id_ = graph.id_;
  return *this;
}

Graph& Graph::operator=(Graph&& graph) {
  this->edge_container_ = move(graph.edge_container_);
  this->nodes_ = move(graph.nodes_);
  this->id_ = move(graph.id_);
  return *this;
}

vector<pair<int, GraphNode>> Graph::getIsolatedNodes(void) {
  vector<pair<int, GraphNode>> isolated_nodes;
  for (const pair<int, GraphNode>& id_and_node : nodes_) {
    if (edge_container_.vertexExist(id_and_node.first)) {
      if (edge_container_.getDegree(id_and_node.first) == 0) {
        pair<int, GraphNode> id_and_node_copy(id_and_node.first,
                                              id_and_node.second);
        isolated_nodes.push_back(id_and_node_copy);
      }
    } else {
      pair<int, GraphNode> id_and_node_copy(id_and_node.first,
                                            id_and_node.second);
      isolated_nodes.push_back(id_and_node_copy);
    }
  }
  return isolated_nodes;
}

vector<int> Graph::getVerticesMissingNodes(void) {
  vector<int> missing;
  vector<int> vertices = edge_container_.getVertices();
  for (int& vertex : vertices) {
    if (nodes_.count(vertex) == 0) {
      missing.push_back(vertex);
    }
  }
  return missing;
}

vector<pair<int, GraphNode>> Graph::getNeighNodes(int vertex) {
  vector<int> neigh_vertices = edge_container_.getNeighVertices(vertex);
  vector<pair<int, GraphNode>> neigh_ids_and_nodes;
  for (int& neigh_vert : neigh_vertices) {
    auto id_and_node = pair<int, GraphNode>(neigh_vert, nodes_[neigh_vert]);
    neigh_ids_and_nodes.push_back(id_and_node);
  }
  return neigh_ids_and_nodes;
}

void Graph::setNode(int vertex, GraphNode graph_node) {
  if (nodes_.count(vertex)) {
    nodes_[vertex] = graph_node;
  } else {
    string errMsg = "Vertex does not exist within graph cannot, reset node";
    throw runtime_error(errMsg);
  }
  calcId_();
}

void Graph::setNode(std::pair<int, GraphNode> id_and_node) {
  setNode(id_and_node.first, id_and_node.second);
}

GraphNode Graph::getNode(int vertex) {
  assert(nodes_.count(vertex));
  return nodes_[vertex];
}

vector<pair<int, GraphNode>> Graph::getNodes(void) {
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

ostream& operator<<(ostream& os, const Graph graph) {
  os << "Graph" << endl;
  for (const pair<int, GraphNode>& id_and_node : graph.nodes_) {
    os << "Node " << id_and_node.first << endl;
    os << id_and_node.second << endl;
  }
  return os;
}

bool cmpVertNodePair(pair<int, GraphNode>& id_and_node1,
                     pair<int, GraphNode>& id_and_node2) {
  string str1_Id = id_and_node1.second.getStringId();
  return str1_Id.compare(id_and_node2.second.getStringId()) < 0;
}
}  // namespace tools
}  // namespace votca
