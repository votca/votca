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

#ifndef _VOTCA_TOOLS_EDGECONTAINER_H
#define _VOTCA_TOOLS_EDGECONTAINER_H

#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <votca/tools/edge.h>

namespace votca {
namespace tools {

/**
 *   \brief EdgeContainer is responsible for operations on groups of edges
 *
 *   Such operations have to deal with finding edges attached to a vertex
 *   finding neighboring vertices etc.
 */
class EdgeContainer {
 protected:
  /// The vertex, the neighboring vertices and the number of edges
  std::unordered_map<int, std::unordered_map<int, int>> adj_list_;

 public:
  /// Constructors can take no arguments a single Edge or a vector of edges
  EdgeContainer() = default;
  ;
  EdgeContainer(Edge edge);
  EdgeContainer(std::vector<Edge> edges);

  /// Get the value of the max degree
  int getMaxDegree() const;
  /// Contains vector of all vertices with degree
  std::vector<int> getVerticesDegree(int degree) const;
  /// Determine the degree of the vertex/number of edges attached
  int getDegree(const int vertex) const;
  /// Determine if a vertex with the specified degree exists
  bool vertexExistWithDegree(int degree) const;

  /// Check if the edge exists returns true or false
  bool edgeExist(const Edge& edge) const;
  /// Check if the vertex exists returns true or false
  bool vertexExist(int vertex) const;
  /// Add an edge to the container
  void addEdge(Edge edge);
  /// Add a lone vertex
  void addVertex(int vertex);
  /// Get all the edges in vector form
  std::vector<Edge> getEdges() const;
  /// Get all the vertices in vector form
  std::vector<int> getVertices() const;
  /// Get the vertices neighboring vert
  std::vector<int> getNeighVertices(int vertex) const;
  /// Get the edges neighboring vert
  std::vector<Edge> getNeighEdges(int vertex) const;
  /// Print output of object
  friend std::ostream& operator<<(std::ostream& os,
                                  const EdgeContainer edge_container);
};
}  // namespace tools
}  // namespace votca

#endif  // _VOTCA_TOOLS_EDGECONTAINER_H
