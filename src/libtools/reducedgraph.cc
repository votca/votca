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

#include <cassert>
#include <votca/tools/edge.h>
#include <votca/tools/reducedgraph.h>

using namespace std;

namespace votca {
namespace tools {

class GraphNode;

bool compareChainWithChains_(const vector<int>& chain,
                             const vector<vector<int>>& chains) {
  bool match = false;
  for (const vector<int>& chain2 : chains) {
    if (chain2.size() == chain.size()) {
      match = true;
      for (size_t index = 0; index < chain.size(); ++index) {
        if (chain2.at(index) != chain.at(index)) {
          match = false;
          break;
        }
      }  // Cycle vertices in each chain
      if (match) return true;
    }  // Chains same size
  }
  return false;
}

void ReducedGraph::init_(vector<ReducedEdge> reduced_edges,
                         unordered_map<int, GraphNode> nodes) {
  vector<Edge> edges;
  nodes_ = nodes;
  for (const ReducedEdge& reduced_edge : reduced_edges) {
    Edge edge(reduced_edge.getEndPoint1(), reduced_edge.getEndPoint2());
    edges.push_back(edge);
    if (expanded_edges_.count(edge)) {
      bool match = compareChainWithChains_(reduced_edge.getChain(),
                                           expanded_edges_[edge]);
      if (!match) {
        expanded_edges_[edge].push_back(reduced_edge.getChain());
      }
    } else {
      expanded_edges_[edge].push_back(reduced_edge.getChain());
    }
  }
  edge_container_ = EdgeContainer(edges);
  calcId_();
}

set<int> getAllVertices_(const std::vector<ReducedEdge>& reduced_edges) {
  set<int> vertices;
  for (const ReducedEdge& reduced_edge : reduced_edges) {
    vector<int> chain = reduced_edge.getChain();
    for (const int vertex : chain) {
      vertices.insert(vertex);
    }
  }
  return vertices;
}

ReducedGraph::ReducedGraph(std::vector<ReducedEdge> reduced_edges) {

  set<int> vertices = getAllVertices_(reduced_edges);
  unordered_map<int, GraphNode> nodes;
  for (const int vertex : vertices) {
    GraphNode gn;
    nodes[vertex] = gn;
  }
  init_(reduced_edges, nodes);
}

ReducedGraph::ReducedGraph(std::vector<ReducedEdge> reduced_edges,
                           unordered_map<int, GraphNode> nodes) {
  set<int> vertices = getAllVertices_(reduced_edges);
  if (nodes.size() < vertices.size()) {
    throw invalid_argument(
        "The number of nodes passed into a reduced graph "
        "must be greater or equivalent to the number of vertices");
  }
  for (const int vertex : vertices) {
    if (nodes.count(vertex) == 0) {
      throw invalid_argument("A vertex is missing its corresponding node.");
    }
  }
  init_(reduced_edges, nodes);
}

Graph ReducedGraph::expandGraph() {
  vector<Edge> all_expanded_edges;
  for (const pair<Edge, vector<vector<int>>> edge_and_vertices :
       expanded_edges_) {
    vector<vector<Edge>> edges_local = expandEdge(edge_and_vertices.first);
    for (vector<Edge>& edges : edges_local) {
      all_expanded_edges.insert(all_expanded_edges.end(), edges.begin(),
                                edges.end());
    }
  }
  return Graph(all_expanded_edges, nodes_);
}

vector<vector<Edge>> ReducedGraph::expandEdge(const Edge& edge) const {
  vector<vector<Edge>> all_edges;
  vector<vector<int>> chains = expanded_edges_.at(edge);
  for (vector<int> vertices : chains) {
    vector<Edge> edges;
    for (size_t index = 0; index < (vertices.size() - 1); ++index) {
      Edge edge_temp(vertices.at(index), vertices.at(index + 1));
      edges.push_back(edge_temp);
    }
    all_edges.push_back(edges);
  }
  return all_edges;
}

set<int> getAllConnectedVertices_(
    const unordered_map<Edge, vector<vector<int>>>& expanded_edges) {
  set<int> all_vertices;

  for (const auto& edge_and_chains : expanded_edges) {
    for (vector<int> chain : edge_and_chains.second) {
      all_vertices.insert(chain.begin(), chain.end());
    }
  }
  return all_vertices;
}
vector<pair<int, GraphNode>> ReducedGraph::getNodes() const {
  vector<int> vertices = edge_container_.getVertices();
  vector<pair<int, GraphNode>> nodes;
  for (const int vertex : vertices) {
    pair<int, GraphNode> id_and_node(vertex, nodes_.at(vertex));
    nodes.push_back(id_and_node);
  }

  set<int> all_connected_vertices = getAllConnectedVertices_(expanded_edges_);
  // Grab the nodes that are not attached to any edges
  for (pair<int, GraphNode> id_and_node : nodes_) {
    if (!all_connected_vertices.count(id_and_node.first)) {
      nodes.push_back(id_and_node);
    }
  }
  return nodes;
}

vector<int> ReducedGraph::getVerticesDegree(int degree) const {
  if (degree == 0) {
    set<int> all_connected_vertices = getAllConnectedVertices_(expanded_edges_);
    vector<int> vertices;
    for (const pair<int, GraphNode> id_and_node : nodes_) {
      if (all_connected_vertices.count(id_and_node.first) == false) {
        vertices.push_back(id_and_node.first);
      }
      return vertices;
    }
  } else {
    return edge_container_.getVerticesDegree(degree);
  }
}

ostream& operator<<(ostream& os, const ReducedGraph graph) {
  os << "Graph" << endl;
  for (const pair<int, GraphNode>& id_and_node : graph.nodes_) {
    os << "Node " << id_and_node.first << endl;
    os << id_and_node.second << endl;
  }

  os << endl;
  os << graph.edge_container_ << endl;

  os << endl;
  os << "Expanded Edge Chains" << endl;

  for (const pair<const Edge, vector<vector<int>>>& edge_and_chains :
       graph.expanded_edges_) {
    for (const vector<int>& chain : edge_and_chains.second) {
      for (const int vertex : chain) {
        os << vertex << " ";
      }
      os << endl;
    }
    os << endl;
  }

  return os;
}

}  // namespace tools
}  // namespace votca
