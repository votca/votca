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

#include <votca/tools/edge.h>

#ifndef _VOTCA_TOOLS_REDUCEDEDGE_H
#define _VOTCA_TOOLS_REDUCEDEDGE_H

namespace votca {
namespace tools {

/**
 * \brief Connects two vertices, also stores the vertices between the endpoints.
 *
 * The edge class stores the ids of two seperate vertices indictating that they
 * are connected (id1,id2). The vertex with the lower id is always at EndPoint1
 * and the vertex with the large id is always placed at end point 2.
 *
 * Unlike, the more basic edge class the reduced edge class allows the user to
 * store a chain of vertices that exist between the two endpoints. For instance
 * an edge of the form:
 *
 * 1 - 2 - 3 - 4
 *
 * Can be stored, in this case EndPoint1 will be vertex 1 and EndPoint2 will be
 * vertex 4. The class name is called a reduced edge because the edge is treated
 * as if it were an edge of only two vertices described as
 *
 * 1 - 4
 *
 **/
class ReducedEdge : public Edge {

 public:
  ReducedEdge() = default;
  /// Creates an edge the smallest integer value will be placed in the id1
  /// spot and the larger in the id2 spot
  ReducedEdge(std::vector<Index> chain);

  ReducedEdge(Index vertex1, Index vertex2)
      : ReducedEdge(std::vector<Index>{vertex1, vertex2}){};

  /// Returns a vector of vertices that constitute the edge
  std::vector<Index> getChain() const { return vertices_; }

  /**
   * \brief Provided a vertex this method will determine if it exists within the
   * reduced edge.
   *
   * Given a reduced edge
   *
   * 1 - 2 - 3 - 4
   *
   * Calling the method with vertex 3 will return true as the vertex exists
   * within the chain stored by the reduced edge.
   *
   * \param[in] Index potential vertex stored in the reduced edge
   * \return bool true if it is in the chain false otherwise
   **/
  bool vertexExistInChain(const int& vertex) const;

  /**
   * \brief Method expands the edge into its parts
   *
   * This method will take an expanded edge e.g.
   *
   * 1 - 2 - 3 - 4 - 5
   *
   * And break it up into its more primitive edges:
   *
   * 1 - 2,
   * 2 - 3,
   * 3 - 4,
   * 4 - 5
   *
   * These edges will then be placed in a vector and returned.
   *
   * \return vector of edges
   **/
  std::vector<Edge> expand() const;

  /// Checks if Edges are equivalent
  bool operator==(const ReducedEdge& edge) const;
  /// Checks if Edges are not equivalent
  bool operator!=(const ReducedEdge& edge) const;

  /**
   * The Reduced Edge follows the same rules as the more basic edge when
   * the EndPoints are differet between the two edges
   *
   * If the vertices are smaller in value
   * Edge ed1(2,3);
   * Edge ed2(1,5);
   * Edge ed3(4,1);
   * priority is given to the smallest vertex
   * assert(ed2<ed1); // will return true
   * assert(ed3<ed1); // will return true
   * assert(ed3<ed1); // will return true
   *
   * In the case that the two edges have the same endpoints:
   *
   * Edge ed4(vector<Index>{ 1, 5, 4});
   * Edge ed5(vector<Index>{ 1, 2, 3, 4});
   *
   * The shorter of the two edges is considered smaller
   * assert(ed4<ed5); // will return true
   *
   * In the case that both edges have the same number of vertices the edge with
   * the first vertex in the chain that is lower in value will be considered
   * smaller.
   *
   * Edge ed6(vector<Index>{ 1, 2, 5, 4});
   * assert(ed5<ed6); // will return true;
   **/

  bool operator<(const ReducedEdge& edge) const;
  bool operator>(const ReducedEdge& edge) const;
  bool operator<=(const ReducedEdge& edge) const;
  bool operator>=(const ReducedEdge& edge) const;

  /// Print the contents of the edge
  friend std::ostream& operator<<(std::ostream& os, const ReducedEdge& edge);
};

}  // namespace tools
}  // namespace votca

/// Define a hasher so we can use it as a key in an unordered_map
namespace std {
template <>
class hash<votca::tools::ReducedEdge> {
 public:
  size_t operator()(const votca::tools::ReducedEdge& ed) const {
    size_t value = 1;
    auto vertices = ed.getChain();
    for (auto vertex : vertices) {
      value += hash<size_t>()(static_cast<size_t>(vertex)) ^ value;
    }
    return value;
  }
};
}  // namespace std
#endif  // _VOTCA_TOOLS_REDUCEDEDGE_H
