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
#ifndef VOTCA_TOOLS_EDGE_H
#define VOTCA_TOOLS_EDGE_H
#pragma once

// Standard includes
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

// Local VOTCA includes
#include "attributes.h"
#include "types.h"

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
class Edge : public Attributes {
 protected:
  std::vector<Index> vertices_;

  //  GraphNode edge_values_;

 public:
  Edge() = default;
  virtual ~Edge() = default;
  /// Creates an edge the smallest integer value will be placed in the id1
  /// spot and the larger in the id2 spot
  Edge(Index ID1, Index ID2);
  /// Given one of the integers in the edge the other will be output
  Index getOtherEndPoint(Index ver) const;
  /// grab the smaller integer
  Index getEndPoint1() const { return vertices_.front(); }
  /// grab the larger integer
  Index getEndPoint2() const { return vertices_.back(); }
  /*
    void addDouble(std::string label, const double& value) {
      edge_values_.addDouble(label, value);
    }

    void addInt(std::string label, const Index value) {
      edge_values_.addInt(label, value);
    }

    void addStr(std::string label, const std::string& value) {
      edge_values_.addStr(label, value);
    }

    double getDouble(const std::string label) const {
      return edge_values_.getDouble(label);
    }

    Index getInt(const std::string label) const {
      return edge_values_.getInt(label);
    }

    std::string getStr(const std::string label) const {
      return edge_values_.getStr(label);
    }

    bool hasInt(const std::string label) const {
      return edge_values_.hasInt(label);
    }*/
  /*
    std::string getEdgePropertiesAsStr() const {
      return edge_values_.getStringId();
    }*/
  /**
   * \brief Checks to see if an edge loops back on itself.
   *
   * If both ends of the edge point to the same vertex than it is considered a
   * loop.
   **/
  bool loop() const { return vertices_.front() == vertices_.back(); }

  /// Determine if the edge contains the Index ID
  bool contains(Index ID) const;
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
const Edge DUMMY_EDGE(std::numeric_limits<Index>::max(),
                      std::numeric_limits<Index>::max());
}  // namespace tools
}  // namespace votca

/// Define a hasher so we can use it as a key in an unordered_map
namespace std {
template <>
class hash<votca::tools::Edge> {
 public:
  size_t operator()(const votca::tools::Edge& ed) const {
    return hash<votca::Index>()(ed.getEndPoint1()) ^
           hash<votca::Index>()(ed.getEndPoint2());
  }
};
}  // namespace std
#endif  // VOTCA_TOOLS_EDGE_H
