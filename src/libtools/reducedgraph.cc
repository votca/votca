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

/******************************************************************************
 * Internal Private Functions
 ******************************************************************************/

/**
 * \brief Function will compare a chain with a vector of chains
 *
 * If the one of the chains in the vector of chains is equivalent will return
 * true. A chain is simply a vector of integers where each integer represents
 * a vertex along an edge. The vertices appear in the chain in the same order
 * as they exist in a graph.
 **/
bool compareChainWithChains_(const vector<Index>& chain,
                             const vector<vector<Index>>& chains) {
  bool match = false;
  for (const vector<Index>& chain2 : chains) {
    if (chain2.size() == chain.size()) {
      match = true;
      for (size_t index = 0; index < chain.size(); ++index) {
        if (chain2.at(index) != chain.at(index)) {
          match = false;
          break;
        }
      }  // Cycle vertices in each chain
      if (match) {
        return true;
      }
    }  // Chains same size
  }
  return false;
}

/**
 * \brief Given a vector of Reduced Edges will find and return junctions if any
 * exist
 *
 * A junction is a vertex in a graph that has 3 or more edges eminating from it.
 *
 * Consider three reduced edges given by
 *
 * Edge   Reduced Edge when Expanded (chain of vertices)
 * 1, 5  :  1 - 3 - 4 - 6 - 5
 * 5, 10 :  5 - 2 - 10
 * 5, 9  :  5 - 8 - 9
 *
 * Graph:
 *
 * 1 - 3 - 4 - 6 - 5 - 2 - 10
 *                 |
 *                 8 - 9
 *
 * Reduced Graph:
 *
 * 1 - 5 - 10
 *     |
 *     9
 *
 * Here vertex 5 would be found to be a junction
 *
 * @param[in,out] - vector of reduced edges
 * @return - set of integers containing junctions
 **/
set<Index> getVertexJunctions_(const vector<ReducedEdge>& reduced_edges) {
  unordered_map<Index, Index> vertex_count;
  for (ReducedEdge reduced_edge : reduced_edges) {
    // if loop, increment value is double and the first index is skipped to
    // prevent over counting of the first number
    Index increment = 1;
    size_t index = 0;
    if (reduced_edge.loop()) {
      ++index;
      increment = 2;
    }
    vector<Index> chain = reduced_edge.getChain();
    for (; index < chain.size(); ++index) {
      if (vertex_count.count(chain.at(index))) {
        vertex_count[chain.at(index)] += increment;
      } else {
        vertex_count[chain.at(index)] = increment;
      }
    }
  }
  set<Index> junctions;
  for (pair<Index, Index> vertex_and_count : vertex_count) {
    if (vertex_and_count.second > 2) {
      junctions.insert(vertex_and_count.first);
    }
  }
  return junctions;
}

/**
 * \breif function adds an edge to an unordered_map and a vector if it is found
 * to not be a loop
 *
 * A loop is an edge that essentially makes a circle. In an edge that is a loop
 * will have both end points equal to the samve vertex. E.g. a reduced edge
 * like this:
 *
 * Edge  Expanded edge (chain of vertices)
 * 2, 2 : 2 - 4 - 3 - 5 - 2
 *
 * Graph
 *
 *  2 - 4
 *  |   |
 *  5 - 3
 *
 * Reduced Graph
 *
 *  - 2
 *  |_|
 *
 * @param[in,out] - vector of edges, an edge is added to the vector if it is not
 * a loop
 * @param[in] - a reduced edge, the edge to be added
 * @param[in,out] - an unordered_map that stores the edge and its chain if it
 * is found to not be a loop
 **/
void addEdgeIfNotLoop_(
    vector<Edge>& edges, const ReducedEdge reduced_edge,
    unordered_map<Edge, vector<vector<Index>>>& expanded_edges) {

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

/**
 * \brief This is a helper function used to help store the vertex chain in a
 * reproducable order
 *
 * Ordering chains in a reproducable manner enable quicker comparison between
 * chains.
 *
 * E.g. Given a chain
 *
 * End Point 1
 * v
 * 1 - 4 - 3 - 5 - 6 - 2 - 1
 *                         ^
 *                      End Point 2
 *
 * This function will ensure that it is reordered such that the next consecutive
 * number in the chain after end point 2 is the lower number, in this case the
 * two numbers that are compared are 4 and 2. The 2 is smaller so the chain is
 * reordered.
 *
 * 1 - 2 - 6 - 5 - 3 - 4 - 1
 *
 * @param[in,out] - vector of integers containing the chain
 **/
void orderChainAfterInitialVertex_(vector<Index>& chain) {
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

bool reorderAndStoreChainIfDoesNotExist_(
    vector<Edge>& edges,
    unordered_map<Edge, vector<vector<Index>>>& expanded_edges,
    vector<Index> chain, Index vertex, size_t& chain_index) {

  Edge edge(vertex, vertex);
  edges.push_back(edge);
  vector<Index> new_chain;
  for (size_t index = 0; index < chain.size(); ++index) {
    if (((chain_index + index) % chain.size()) == 0) {
      ++chain_index;
    }
    Index new_chain_index = (chain_index + index) % chain.size();
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

set<Index> getAllVertices_(const std::vector<ReducedEdge>& reduced_edges) {
  set<Index> vertices;
  for (const ReducedEdge& reduced_edge : reduced_edges) {
    vector<Index> chain = reduced_edge.getChain();
    for (const Index vertex : chain) {
      vertices.insert(vertex);
    }
  }
  return vertices;
}

set<Index> getAllConnectedVertices_(
    const unordered_map<Edge, vector<vector<Index>>>& expanded_edges) {
  set<Index> all_vertices;

  for (const auto& edge_and_chains : expanded_edges) {
    for (vector<Index> chain : edge_and_chains.second) {
      all_vertices.insert(chain.begin(), chain.end());
    }
  }
  return all_vertices;
}
/******************************************************************************
 * Private Class Methods
 ******************************************************************************/

void ReducedGraph::init_(vector<ReducedEdge> reduced_edges,
                         unordered_map<Index, GraphNode> nodes) {
  vector<Edge> edges;
  nodes_ = nodes;

  junctions_ = getVertexJunctions_(reduced_edges);

  for (const ReducedEdge& reduced_edge : reduced_edges) {

    if (reduced_edge.loop() &&
        junctions_.count(reduced_edge.getEndPoint1()) == 0) {
      vector<Index> chain = reduced_edge.getChain();
      size_t chain_index = 0;
      bool edge_added = false;
      for (Index vertex : chain) {
        if (junctions_.count(vertex)) {
          edge_added = reorderAndStoreChainIfDoesNotExist_(
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

/******************************************************************************
 * Public Class Methods
 ******************************************************************************/

ReducedGraph::ReducedGraph(std::vector<ReducedEdge> reduced_edges) {

  set<Index> vertices = getAllVertices_(reduced_edges);
  unordered_map<Index, GraphNode> nodes;
  for (const Index vertex : vertices) {
    GraphNode gn;
    nodes[vertex] = gn;
  }
  init_(reduced_edges, nodes);
}

ReducedGraph::ReducedGraph(std::vector<ReducedEdge> reduced_edges,
                           unordered_map<Index, GraphNode> nodes) {

  set<Index> vertices = getAllVertices_(reduced_edges);
  if (nodes.size() < vertices.size()) {
    throw invalid_argument(
        "The number of nodes passed into a reduced graph "
        "must be greater or equivalent to the number of vertices");
  }
  for (const Index vertex : vertices) {
    if (nodes.count(vertex) == 0) {
      throw invalid_argument("A vertex is missing its corresponding node.");
    }
  }
  init_(reduced_edges, nodes);

  // Add all the nodes that are isolated
  for (pair<const Index, GraphNode>& id_and_node : nodes) {
    if (vertices.count(id_and_node.first) == 0) {
      edge_container_.addVertex(id_and_node.first);
    }
  }
}

Graph ReducedGraph::expandGraph() const {
  vector<Edge> all_expanded_edges;
  for (const pair<Edge, vector<vector<Index>>> edge_and_vertices :
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
  vector<vector<Index>> chains = expanded_edges_.at(edge);
  for (vector<Index> vertices : chains) {
    vector<Edge> edges;
    for (size_t index = 0; index < (vertices.size() - 1); ++index) {
      Edge edge_temp(vertices.at(index), vertices.at(index + 1));
      edges.push_back(edge_temp);
    }
    all_edges.push_back(edges);
  }
  return all_edges;
}

vector<pair<Index, GraphNode>> ReducedGraph::getNodes() const {
  vector<Index> vertices = edge_container_.getVertices();
  vector<pair<Index, GraphNode>> nodes;
  for (const Index vertex : vertices) {
    pair<Index, GraphNode> id_and_node(vertex, nodes_.at(vertex));
    nodes.push_back(id_and_node);
  }

  set<Index> all_connected_vertices = getAllConnectedVertices_(expanded_edges_);
  // Grab the nodes that are not attached to any edges
  for (pair<Index, GraphNode> id_and_node : nodes_) {
    if (!all_connected_vertices.count(id_and_node.first)) {
      nodes.push_back(id_and_node);
    }
  }
  return nodes;
}

vector<Index> ReducedGraph::getVerticesDegree(Index degree) const {
  if (degree == 0) {
    set<Index> all_connected_vertices =
        getAllConnectedVertices_(expanded_edges_);
    vector<Index> vertices;
    for (const pair<Index, GraphNode> id_and_node : nodes_) {
      if (all_connected_vertices.count(id_and_node.first) == false) {
        vertices.push_back(id_and_node.first);
      }
      return vertices;
    }
  }
  return edge_container_.getVerticesDegree(degree);
}

ostream& operator<<(ostream& os, const ReducedGraph graph) {
  os << "Graph" << endl;
  for (const pair<Index, GraphNode>& id_and_node : graph.nodes_) {
    os << "Node " << id_and_node.first << endl;
    os << id_and_node.second << endl;
  }

  os << endl;
  os << graph.edge_container_ << endl;

  os << endl;
  os << "Expanded Edge Chains" << endl;

  for (const pair<const Edge, vector<vector<Index>>>& edge_and_chains :
       graph.expanded_edges_) {
    for (const vector<Index>& chain : edge_and_chains.second) {
      for (const Index vertex : chain) {
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
