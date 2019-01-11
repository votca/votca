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

#include <votca/tools/edge.h>

#ifndef _VOTCA_TOOLS_REDUCEDEDGE_H
#define _VOTCA_TOOLS_REDUCEDEDGE_H

namespace votca {
namespace tools {

/**
 * \brief Connects to vertices
 *
 * The edge class stores the ids of two seperate vertices indictating that they
 * are connected (id1,id2). Unlike a pair the vertex with the lower value is
 * always placed in id1, this allows us to reduce ambiguity when dealing with
 * a link.
 */
class ReducedEdge : public Edge{

 public:
  ReducedEdge() {}
  ~ReducedEdge() {}
  /// Creates an edge the smallest integer value will be placed in the id1
  /// spot and the larger in the id2 spot
  ReducedEdge(std::vector<int> chain);
  ReducedEdge(int vert1, int vert2) : ReducedEdge(std::vector<int>{vert1,vert2}) {};

  std::vector<int> getChain() const { return vertices_;}

  std::vector<Edge> expand() const;
  /// Checks if Edges are equivalent
  bool operator==(const ReducedEdge ed) const;
  /// Checks if Edges are not equivalent
  bool operator!=(const ReducedEdge ed) const;

  /// If the vertices are smaller in value
  /// Edge ed1(2,3);
  /// Edge ed2(1,5);
  /// Edge ed3(4,1);
  /// priority is given to the smallest vertex
  /// assert(ed2<ed1); // will return true
  /// assert(ed3<ed1); // will return true
  /// assert(ed3<ed1); // will return true
  bool operator<(const ReducedEdge ed) const;
  bool operator>(const ReducedEdge ed) const;
  bool operator<=(const ReducedEdge ed) const;
  bool operator>=(const ReducedEdge ed) const;

  /// Print the contents of the edge
  friend std::ostream& operator<<(std::ostream& os, const ReducedEdge ed);
};

// Value used as a dummy object
const ReducedEdge DUMMY_REDUCEDEDGE(
    std::vector<int>{std::numeric_limits<int>::max(),std::numeric_limits<int>::max()});
}
}

/// Define a hasher so we can use it as a key in an unordered_map
namespace std {
template <>
class hash<votca::tools::ReducedEdge> {
 public:
  size_t operator()(const votca::tools::ReducedEdge& ed) const {
    size_t value = 1;
    auto vertices = ed.getChain();
    for(auto vertex : vertices ){
      value+=hash<size_t>()(static_cast<size_t>(vertex))^value;
    }
    return value;
  }
};
}
#endif  // _VOTCA_TOOLS_REDUCEDEDGE_H
