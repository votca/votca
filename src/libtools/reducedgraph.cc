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

set<int> getVertexJunctions_(const vector<ReducedEdge>& reduced_edges) {
  unordered_map<int, int> vertex_count;
  for (ReducedEdge reduced_edge : reduced_edges) {
    // if loop, increment value is double and the first index is skipped to
    // prevent over counting of the first number
    int increment = 1;
    size_t index = 0;
    if (reduced_edge.loop()) {
      ++index;
      increment = 2;
    }
    vector<int> chain = reduced_edge.getChain();
    for (; index < chain.size(); ++index) {
      if (vertex_count.count(chain.at(index))) {
        vertex_count[chain.at(index)] += increment;
      } else {
        vertex_count[chain.at(index)] = increment;
      }
    }
  }
  set<int> junctions;
  for (pair<int, int> vertex_and_count : vertex_count) {
    if (vertex_and_count.second > 2) {
      junctions.insert(vertex_and_count.first);
    }
  }
  return junctions;
}

void addEdgeIfNotLoop_(
    vector<Edge>& edges, const ReducedEdge reduced_edge,
    unordered_map<Edge, vector<vector<int>>>& expanded_edges) {

  Edge edge(reduced_edge.getEndPoint1(), reduced_edge.getEndPoint2());
  if (expanded_edges.count(edge)) {
    bool match =
        compareChainWithChains_(reduced_edge.getChain(), expanded_edges[edge]);
    if (!match) {
      expanded_edges[edge].push_back(reduced_edge.getChain());
    }
  } else {
    expanded_edges[edge].push_back(reduced_edge.getChain());
  }
  edges.push_back(edge);
}

void orderChainAfterInitialVertex_(vector<int>& chain) {
  size_t ignore_first_and_last_vertex = 2;
  size_t total_number_to_parse =
      (chain.size() - ignore_first_and_last_vertex) / 2;
  bool reverse_vector = false;
  for (size_t count = 0; count < (total_number_to_parse); ++count) {
    if (chain.at(chain.size() - count - 1) < chain.at(count + 1)) {
      reverse_vector = true;
      break;
    }
  }

  if (reverse_vector) {
    reverse(chain.begin(), chain.end());
  }
}

bool reordereAndStoreChainIfDoesNotExist_(
    vector<Edge>& edges,
    unordered_map<Edge, vector<vector<int>>>& expanded_edges, vector<int> chain,
    int vertex, size_t& chain_index) {

  Edge edge(vertex, vertex);
  edges.push_back(edge);
  vector<int> new_chain;
  for (size_t index = 0; index < chain.size(); ++index) {
    if (((chain_index + index) % chain.size()) == 0) {
      ++chain_index;
    }
    int new_chain_index = (chain_index + index) % chain.size();
    new_chain.push_back(chain.at(new_chain_index));
  }
  // Ensure that the new_chain is sorted so after the first vertex they are
  // ordered from smallest to largest
  orderChainAfterInitialVertex_(new_chain);
  bool match = compareChainWithChains_(new_chain, expanded_edges[edge]);
  if (!match) {
    expanded_edges[edge].push_back(new_chain);
    return true;
  }
  return false;
}

void ReducedGraph::init_(vector<ReducedEdge> reduced_edges,
                         unordered_map<int, GraphNode> nodes) {
  vector<Edge> edges;
  nodes_ = nodes;

  junctions_ = getVertexJunctions_(reduced_edges);

  for (const ReducedEdge& reduced_edge : reduced_edges) {

    if (reduced_edge.loop() &&
        junctions_.count(reduced_edge.getEndPoint1()) == 0) {
      vector<int> chain = reduced_edge.getChain();
      size_t chain_index = 0;
      bool edge_added = false;
      for (int vertex : chain) {
        if (junctions_.count(vertex)) {
          edge_added = reordereAndStoreChainIfDoesNotExist_(
              edges, expanded_edges_, chain, vertex, chain_index);
          break;
        }
        ++chain_index;
      }
      if (!edge_added) {
        addEdgeIfNotLoop_(edges, reduced_edge, expanded_edges_);
      }
    } else {
      addEdgeIfNotLoop_(edges, reduced_edge, expanded_edges_);
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

  // Add all the nodes that are isolated
  for (pair<const int, GraphNode>& id_and_node : nodes) {
    if (vertices.count(id_and_node.first) == 0) {
      edge_container_.addVertex(id_and_node.first);
    }
  }
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
/*
bool ReducedGraph::edgeExist(const Edge& edge) const {
  return edge_container_.edgeExist(edge);
}*/

vector<Edge> ReducedGraph::getEdges() {
  /*  vector<Edge> edges_unfiltered = edge_container_.getEdges();
    vector<Edge> edges;
    for(const Edge edge : edges_unfiltered){
     if(edge.loop()){
      // Find which vertex if any is a junction if it is return the edge so that
      // it points to the junction
      vector<vector<int>> chains = expanded_edges_.at(edge);
      // If it is a loop there should only be a single junction at most
      assert(chains.size()==1);
      for(int vertex : chains.at(0)){
        if(junctions_.count(vertex)){
          Edge edge(vertex,vertex);
          edges.push_back(edge);
          break;
        }
      }
     }else{
      edges.push_back(edge);
     }
     }
    return edges;*/
  return edge_container_.getEdges();
}

vector<int> ReducedGraph::getVertices() const {
  // Get all the vertices that are in the reduced graph
  vector<int> vertices = edge_container_.getVertices();
  return vertices;
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
