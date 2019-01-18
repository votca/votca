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

#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <votca/tools/edgecontainer.h>
#include <votca/tools/graphnode.h>

#ifndef _VOTCA_TOOLS_GRAPH_H
#define _VOTCA_TOOLS_GRAPH_H

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

  /// Parameter description
  /// @param int - is the index of the graph nodes / vertex ids
  /// @param GraphNode - this is the node object at each vertex and contains
  /// all the informatino that is relevant to that object
  std::unordered_map<int, GraphNode> nodes_;

  /// This is the id of the graph to graphs that contain the same content
  /// are considered equal
  std::string id_;

 protected:
  /// Calculate the id of the graph
  void calcId_();

 public:
  Graph() : id_(""){};
  virtual ~Graph(){};
  /// Constructor
  /// @param edgs - vector of edges where each edge is composed of two
  /// ints (vertex ids) describing a link between the vertices
  /// @param nodes - unordered_map where the key is the vertex id and the
  /// target is the graph node
  Graph(std::vector<Edge> edgs, std::unordered_map<int, GraphNode> nodes)
      : edge_container_(edgs), nodes_(nodes) {
    calcId_();
  }

  /// Equivalence and non equivalence operators work by determine if the
  /// contents of each graph node in each of the graphs are the same.
  bool operator!=(const Graph& graph) const;
  bool operator==(const Graph& graph) const;

  /// Find all the vertices that are isolated (not connected to any other
  /// vertex) and return them in a vector with their corresponding graph node.
  std::vector<std::pair<int, GraphNode>> getIsolatedNodes(void) const;

  /// Functions determines which vertices do not have a graph node associated
  /// with them and return their ids in a vector.
  std::vector<int> getVerticesMissingNodes(void);

  /// Returns a vector of the vertices and their graph nodes that are directly
  /// connected to the vertex 'vert'
  std::vector<std::pair<int, GraphNode>> getNeighNodes(int vertex) const;

  /// set the Node associated with vertex 'vert'
  void setNode(int vertex, GraphNode& graph_node);
  void setNode(std::pair<int, GraphNode>& id_and_node);

  /// Gets all vertices with degree of 3 or greater
  std::vector<int> getJunctions() const;

  /// Return a copy of the graph node at vertex 'vert'
  GraphNode getNode(const int vertex) const;

  /// Return all the vertices and their graph nodes that are within the graph
  virtual std::vector<std::pair<int, GraphNode>> getNodes(void) const;

  /// Returns all the vertices of the graph connected to vertex `vert` through
  /// an edge.
  std::vector<int> getNeighVertices(int vertex) const {
    return edge_container_.getNeighVertices(vertex);
  }

  /// Returns the id of graph
  std::string getId() const { return id_; }

  /// Returns all the edges in the graph
  std::vector<Edge> getEdges() { return edge_container_.getEdges(); }

  /// Returns all the edges in the graph connected to vertex `vertex`
  std::vector<Edge> getNeighEdges(int vertex) const {
    return edge_container_.getNeighEdges(vertex);
  }

  /// Returns all the vertices in the graph
  std::vector<int> getVertices() const { return edge_container_.getVertices(); }

  /**
   * \brief Finds the max degree of a vertex in the graph.
   *
   * Will look at each of the vertices and find the vertex with the most edges
   * connected to it. It will count the number of edges this corresponds to the
   * maximum degree of the graph.
   **/
  int getMaxDegree() const { return edge_container_.getMaxDegree(); }

  /// Calcualtes the degree, or number of edges connected to vertex `vertex`
  int getDegree(int vertex) const;

  /// Returns all the vertices with degree specified by `degree`
  virtual std::vector<int> getVerticesDegree(int degree) const;

  /// Determines if a vertex exists within the graph
  bool vertexExist(int vertex) const;

  /// Determines if an edge is stored in the graph
  bool edgeExist(const Edge& edge) const {
    return edge_container_.edgeExist(edge);
  }

  /// Copies nodes from one graph to another. This should only be used in cases
  /// where the graph does not contain nodes before the copy.
  void copyNodes(Graph& graph);

  friend std::ostream& operator<<(std::ostream& os, const Graph graph);
};

/**
 * \brief Compare function pair<int,GraphNode> object
 *
 * This function is meant to be used with the stl sort algorithm. It will sort
 * a vector of pairs containing the vertex ids and the graphnodes. Only the
 * contetns of the graph node object are used to determine precidence e.g.
 *
 * pair<int,GraphNode> pr_grn1{ 1, gn };
 * pair<int,GraphNode> pr_grn2{ 2, gn2 };
 *
 * vector<pair<int,GraphNode> > vec_pr_gn = { pr_grn1, pr_grn2 , ... etc };
 *
 * sort(vec_pr_gn.begin(),vec_pr_gn.end(),cmpVertNodePair);
 */
bool cmpVertNodePair(std::pair<int, GraphNode>& id_and_node1,
                     std::pair<int, GraphNode>& id_and_node2);
}  // namespace tools
}  // namespace votca
#endif  // _VOTCA_TOOLS_GRAPH_H
