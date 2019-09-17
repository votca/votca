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
#include <limits>
#include <utility>
#include <vector>

#ifndef _VOTCA_TOOLS_EDGE_H
#define _VOTCA_TOOLS_EDGE_H

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
class Edge {
 protected:
  std::vector<int> vertices_;

 public:
  Edge() {}
  virtual ~Edge() {}
  /// Creates an edge the smallest integer value will be placed in the id1
  /// spot and the larger in the id2 spot
  Edge(int ID1, int ID2);
  /// Given one of the integers in the edge the other will be output
  int getOtherEndPoint(int ver) const;
  /// grab the smaller integer
  int getEndPoint1() const { return vertices_.front(); }
  /// grab the larger integer
  int getEndPoint2() const { return vertices_.back(); }

  /**
   * \brief Checks to see if an edge loops back on itself.
   *
   * If both ends of the edge point to the same vertex than it is considered a
   * loop.
   **/
  bool loop() const { return vertices_.front() == vertices_.back(); }

  /// Determine if the edge contains the int ID
  bool contains(int ID) const;
  /// Checks if Edges are equivalent
  virtual bool operator==(const Edge& ed) const;
  /// Checks if Edges are not equivalent
  virtual bool operator!=(const Edge& ed) const;

  /// If the vertices are smaller in value
  /// Edge ed1(2,3);
  /// Edge ed2(1,5);
  /// Edge ed3(4,1);
  /// priority is given to the smallest vertex
  /// assert(ed2<ed1); // will return true
  /// assert(ed3<ed1); // will return true
  /// assert(ed3<ed1); // will return true
  virtual bool operator<(const Edge& ed) const;
  virtual bool operator>(const Edge& ed) const;
  virtual bool operator<=(const Edge& ed) const;
  virtual bool operator>=(const Edge& ed) const;

  /// Print the contents of the edge
  friend std::ostream& operator<<(std::ostream& os, const Edge& ed);
};

// Value used as a dummy object
const Edge DUMMY_EDGE(std::numeric_limits<int>::max(),
                      std::numeric_limits<int>::max());
}  // namespace tools
}  // namespace votca

/// Define a hasher so we can use it as a key in an unordered_map
namespace std {
template <>
class hash<votca::tools::Edge> {
 public:
  size_t operator()(const votca::tools::Edge& ed) const {
    return hash<int>()(ed.getEndPoint1()) ^ hash<int>()(ed.getEndPoint2());
  }
};
}  // namespace std
#endif  // _VOTCA_TOOLS_EDGE_H
