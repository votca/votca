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

#include <unordered_map>
#include <string>
#include <iostream>
#include <utility>
#include <votca/tools/name.h>

#ifndef _VOTCA_TOOLS_GRAPHNODE_H
#define _VOTCA_TOOLS_GRAPHNODE_H

namespace votca {
namespace tools {

class GraphDistVisitor;
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
class GraphNode {
 private:
  std::string str_id_{""};
  std::unordered_map<std::string, int> int_vals_;
  std::unordered_map<std::string, double> double_vals_;
  std::unordered_map<std::string, std::string> str_vals_;
  void initStringId_();
 public:
  GraphNode() {};
  /// Constructor
  GraphNode(const std::unordered_map<std::string, int> int_vals,
            const std::unordered_map<std::string, double> double_vals,
            const std::unordered_map<std::string, std::string> str_vals);
  /// Basic setters
  void setInt(const std::unordered_map<std::string, int> int_vals);
  void setDouble(const std::unordered_map<std::string, double> double_vals);
  void setStr(const std::unordered_map<std::string, std::string> str_vals);
  /// Get the string id unique to the contents
  std::string getStringId() const { return str_id_; }
  bool operator==(const GraphNode gn) const;
  bool operator!=(const GraphNode gn) const;
  // Allow visitor to directly access members of the node
  friend GraphDistVisitor;
};

}
}
#endif  // _VOTCA_TOOLS_GRAPHNODE_H
