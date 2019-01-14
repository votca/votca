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
  for (auto reduced_edge : reduced_edges) {
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

vector<vector<Edge>> ReducedGraph::expandEdge(Edge edge) {
  assert(expanded_edges_.count(edge));
  vector<vector<Edge>> all_edges;
  vector<vector<int>> chains = expanded_edges_[edge];
  for (vector<int> vertices : chains) {
    vector<Edge> edges;
    for (size_t index = 0; index < (vertices.size() - 1); ++index) {
      Edge edge(vertices.at(index), vertices.at(index + 1));
      edges.push_back(edge);
    }
    all_edges.push_back(edges);
  }
  return all_edges;
}

vector<pair<int, GraphNode>> ReducedGraph::getNodes() {
  vector<int> vertices = edge_container_.getVertices();
  vector<pair<int, GraphNode>> nodes;
  for (const int vertex : vertices) {
    pair<int, GraphNode> id_and_node(vertex, nodes_[vertex]);
    nodes.push_back(id_and_node);
  }
  return nodes;
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
