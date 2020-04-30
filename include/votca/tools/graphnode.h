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
#ifndef VOTCA_TOOLS_GRAPHNODE_H
#define VOTCA_TOOLS_GRAPHNODE_H

#include "attributes.h"
#include "types.h"

#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>

namespace votca {
namespace tools {

// class GraphDistVisitor;

// enum class StringIdType { full, brief };

/**
 * \brief A graph node that will take a variety of different values
 *
 * The purpose of this object is to allow a level of flexibility when using
 * graph type algorithms. The object uses its contents to create a string
 * that is unique to the contents. If two nodes with the same contents are
 * created they will be considered to be equivalent.
 *
 * NOTE: It may be of interest to take a look at the the Boost property map
 * class, which was designed for a similar purpose.
 */
class GraphNode : public Attributes {
 private:
  // std::string str_id_{""};
  // std::string str_id_no_label_{""};
  // std::unordered_map<std::string, Index> int_vals_;
  // std::unordered_map<std::string, double> double_vals_;
  // std::unordered_map<std::string, std::string> str_vals_;
  // void initStringId_();

 public:
  GraphNode() = default;

  /// Constructor
  /// Each map corresponds to a different content the graph node can contain.
  /* GraphNode(const std::unordered_map<std::string, Index> int_vals,
             const std::unordered_map<std::string, double> double_vals,
             const std::unordered_map<std::string, std::string> str_vals);*/

  template <typename T>
  GraphNode(const std::unordered_map<std::string, T> values)
      : Attributes(values){};

  template <typename T, typename U>
  GraphNode(const std::unordered_map<std::string, T> values1,
            const std::unordered_map<std::string, U> values2)
      : Attributes(values1, values2){};

  template <typename T, typename U, typename V>
  GraphNode(const std::unordered_map<std::string, T> values1,
            const std::unordered_map<std::string, U> values2,
            const std::unordered_map<std::string, V> values3)
      : Attributes(values1, values2, values3){};

  /// Basic setters
  //  void setInt(const std::unordered_map<std::string, Index> int_vals);
  // void setDouble(const std::unordered_map<std::string, double> double_vals);
  // void setStr(const std::unordered_map<std::string, std::string> str_vals);

  // void resetInt(std::string label, const Index& value);
  // void resetDouble(std::string label, const double& value);
  // void resetStr(std::string label, const std::string& value);

  // void addInt(std::string label, const Index value);
  // void addDouble(std::string label, const double value);
  // void addStr(std::string label, const std::string& value);

  // bool hasInt(const std::string& label) const;

  /// Basic getters
  // Index getInt(const std::string& str) const;
  // double getDouble(const std::string& str) const;
  // std::string getStr(const std::string& str) const;

  /// Get the string id unique to the contents of the graph node, by default
  /// will return the full string id
  //  std::string getStringId(StringIdType type = StringIdType::full) const {
  //    if (type == StringIdType::full) return str_id_;
  //    return str_id_no_label_;
  //  }
  // bool operator==(const GraphNode gn) const;
  // bool operator!=(const GraphNode gn) const;

  // Allow visitor to directly access members of the node
  // friend GraphDistVisitor;
  friend std::ostream& operator<<(std::ostream& os, const GraphNode gn);
};

/**
 * \brief Comparison function to be used with stl sort algorithm
 *
 * Given a vector of graph node objects this function can be used to sort them
 * they are sorted based on the contents of the graphnode which are summarized
 * in the graph nodes string id. E.g.
 *
 * vector<GraphNode> vec_gn = { gn1, gn2,...etc };
 *
 * sort(vec_gn.begin(),vec_gn.end(),cmpNode);
 */
bool cmpNode(GraphNode gn1, GraphNode gn2);

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_GRAPHNODE_H
