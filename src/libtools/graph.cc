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
#include <string>
#include <votca/tools/graph.h>

using namespace std;

namespace votca {
namespace tools {

class GraphNode;

void Graph::updateIds_(Graph& g) {
  if (!this->id_set_) {
    this->calcId_();
    this->id_set_ = true;
  }
  if (!g.id_set_) {
    g.calcId_();
    g.id_set_ = true;
  }
}

bool Graph::operator!=(Graph& g) {
  updateIds_(g);
  return id_.compare(g.id_);
}

bool Graph::operator==(Graph& g) { return !(*(this) != g); }

Graph& Graph::operator=(const Graph& g) {
  this->adj_list_ = g.adj_list_;
  for (auto pr : g.nodes_) {
    this->nodes_[pr.first] = pr.second;
  }
  this->id_ = g.id_;
  return *this;
}

vector<pair<int, GraphNode>> Graph::getIsolatedNodes(void) {
  vector<pair<int, GraphNode>> iso_nodes;
  for (auto node : nodes_) {
    if (adj_list_.count(node.first)) {
      if (adj_list_[node.first].size() == 0) {
        pair<int, GraphNode> pr(node.first, node.second);
        iso_nodes.push_back(pr);
      }
    } else {
      pair<int, GraphNode> pr(node.first, node.second);
      iso_nodes.push_back(pr);
    }
  }
  return iso_nodes;
}

vector<int> Graph::getVerticesMissingNodes(void) {
  vector<int> missing;
  for (auto pr_v : adj_list_) {
    if (nodes_.count(pr_v.first) == 0) {
      missing.push_back(pr_v.first);
    }
  }
  return missing;
}

GraphNode& Graph::Node(int vert) { return nodes_[vert]; }

GraphNode Graph::getNode(int vert) { return nodes_[vert]; }

vector<pair<int, GraphNode>> Graph::getNodes(void) {
  vector<pair<int, GraphNode>> vec_nodes;
  for (auto pr_node : nodes_) {
    vec_nodes.push_back(pr_node);
  }
  return vec_nodes;
}

void Graph::calcId_() {
  auto nodes = getNodes();
  sort(nodes.begin(), nodes.end(), cmpVertNodePair);
  string struct_Id_temp = "";
  for (auto nd_pr : nodes) {
    struct_Id_temp.append(nd_pr.second.getStringId());
  }
  id_ = struct_Id_temp;
  return;
}

ostream& operator<<(ostream& os, const Graph g) {
  os << "Graph" << endl;
  for (auto p_gn : g.nodes_) {
    os << "Node " << p_gn.first << endl;
    os << p_gn.second << endl;
  }
  return os;
}

string Graph::getId(void) {
  if (!this->id_set_) {
    this->calcId_();
    this->id_set_ = true;
  }
  return id_;
}

bool cmpVertNodePair(pair<int, GraphNode> gn1_pr, pair<int, GraphNode> gn2_pr) {
  string str1_Id = gn1_pr.second.getStringId();
  return str1_Id.compare(gn2_pr.second.getStringId()) < 0;
}
}
}
