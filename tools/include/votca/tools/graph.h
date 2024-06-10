/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

#ifndef VOTCA_TOOLS_GRAPH_H
#define VOTCA_TOOLS_GRAPH_H

// Standard includes
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Local VOTCA includes
#include "edgecontainer.h"
#include "graphnode.h"

namespace votca {
namespace tools {

/**
 * \brief A graph object that contains the graph nodes and the edges describing
 * the bonds between nodes.
 *
 */
class GraphNode;

class Graph {
 protected:
  EdgeContainer edge_container_;

  std::unordered_map<Index, GraphNode> nodes_;

  /// This is the id of the graph to graphs that contain the same content
  /// are considered equal
  std::string id_;

 protected:
  /// Calculate the id of the graph
  void calcId_();

 public:
  Graph() : id_("") {};
  virtual ~Graph() = default;
  /// Constructor
  /// @param edges - vector of edges where each edge is composed of two
  /// s (vertex ids) describing a link between the vertices
  /// @param nodes - unordered_map where the key is the vertex id and the
  /// target is the graph node
  Graph(std::vector<Edge> edges, std::unordered_map<Index, GraphNode> nodes);

  /// Equivalence and non equivalence operators work by determine if the
  /// contents of each graph node in each of the graphs are the same.
  bool operator!=(const Graph& graph) const;
  bool operator==(const Graph& graph) const;

  /// Find all the vertices that are isolated (not connected to any other
  /// vertex) and return them in a vector with their corresponding graph node.
  std::vector<std::pair<Index, GraphNode>> getIsolatedNodes(void) const;

  /// Returns a vector of the vertices and their graph nodes that are directly
  /// connected to the vertex 'vert'
  std::vector<std::pair<Index, GraphNode>> getNeighNodes(Index vertex) const;

  /// set the Node associated with vertex 'vert'
  void setNode(Index vertex, GraphNode& graph_node);
  void setNode(std::pair<Index, GraphNode>& id_and_node);

  /// Gets all vertices with degree of 3 or greater
  std::vector<Index> getJunctions() const;

  /// Return a copy of the graph node at vertex 'vert'
  GraphNode getNode(const Index vertex) const;

  /// Return all the vertices and their graph nodes that are within the graph
  virtual std::vector<std::pair<Index, GraphNode>> getNodes() const;

  /// Returns all the vertices of the graph connected to vertex `vert` through
  /// an edge.
  std::vector<Index> getNeighVertices(Index vertex) const {
    return edge_container_.getNeighVertices(vertex);
  }

  /// Returns the id of graph
  std::string getId() const { return id_; }

  /// Returns all the edges in the graph
  virtual std::vector<Edge> getEdges() { return edge_container_.getEdges(); }

  /// Returns all the edges in the graph connected to vertex `vertex`
  std::vector<Edge> getNeighEdges(Index vertex) const {
    return edge_container_.getNeighEdges(vertex);
  }

  /// Returns all the vertices in the graph
  std::vector<Index> getVertices() const;

  /**
   * \brief Finds the max degree of a vertex in the graph.
   *
   * Will look at each of the vertices and find the vertex with the most edges
   * connected to it. It will count the number of edges this corresponds to the
   * maximum degree of the graph.
   **/
  Index getMaxDegree() const { return edge_container_.getMaxDegree(); }

  /// Calcualtes the degree, or number of edges connected to vertex `vertex`
  Index getDegree(Index vertex) const;

  /// Returns all the vertices with degree specified by `degree`
  virtual std::vector<Index> getVerticesDegree(Index degree) const;

  /// Determines if a vertex exists within the graph
  bool vertexExist(Index vertex) const;

  /// Determines if an edge is stored in the graph
  virtual bool edgeExist(const Edge& edge) const {
    return edge_container_.edgeExist(edge);
  }

  /// Remove contents of all nodes
  void clearNodes();
  /// Copies nodes from one graph to another. This should only be used in cases
  /// where the graph does not contain nodes before the copy.
  void copyNodes(Graph& graph);

  friend std::ostream& operator<<(std::ostream& os, const Graph graph);
};

/**
 * \brief Compare function pair<Index ,GraphNode> object
 *
 * This function is meant to be used with the stl sort algorithm. It will sort
 * a vector of pairs containing the vertex ids and the graphnodes. Only the
 * contetns of the graph node object are used to determine precidence e.g.
 *
 * pair<Index ,GraphNode> pr_grn1{ 1, gn };
 * pair<Index ,GraphNode> pr_grn2{ 2, gn2 };
 *
 * vector<pair<Index ,GraphNode> > vec_pr_gn = { pr_grn1, pr_grn2 , ... etc };
 *
 * sort(vec_pr_gn.begin(),vec_pr_gn.end(),cmpVertNodePair);
 */
bool cmpVertNodePair(const std::pair<Index, GraphNode>& id_and_node1,
                     const std::pair<Index, GraphNode>& id_and_node2);
}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_GRAPH_H
