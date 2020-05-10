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

#ifndef VOTCA_TOOLS_CONTENTLABEL_H
#define VOTCA_TOOLS_CONTENTLABEL_H
#pragma once

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/any.hpp>
namespace votca {
namespace tools {

class GraphNode;
class Branch;

enum class LabelType { verbose, concise };

/**
 * @brief Content Label is meant to be a unique identifier that contains
 * the contents of some object
 *
 * The content label encodes all of the structural information about an
 * object into a single string.  The idea behind the content label is to
 * generate a label that is unique the structure of a given graph object,
 * such that the objects structure could be simply rebuilt with the content
 * label alone.  Store the contents as a vector of strings that appear in
 * the correct sequence
 */
class ContentLabel {
 public:
  static constexpr int arr_len = 4;
  typedef std::array<std::string, arr_len> KeyValType;

 private:
  // index 0 - opening type, 1 key, 2 - val, 3 - closing type
  std::vector<KeyValType> labels_;
  size_t hash_ = 0;
  size_t label_char_len_ = 0;

  void initLabels_(std::unordered_map<std::string, boost::any> values);
  bool containsBranch_() const;

 public:
  ContentLabel() = default;
  ContentLabel(std::unordered_map<std::string, boost::any> values);
  void add(GraphNode gn);
  void add(Branch br);
  void append(ContentLabel);

  bool isEmpty() const noexcept { return labels_.size() == 0; }
  void clear();

  void makeBranch();

  std::string get(const LabelType& label_type = LabelType::verbose) const;

  bool operator!=(const ContentLabel& label) const;
  bool operator==(const ContentLabel& label) const;
  bool operator<(const ContentLabel& label) const;
  bool operator<=(const ContentLabel& label) const;
  bool operator>(const ContentLabel& label) const;
  bool operator>=(const ContentLabel& label) const;
};

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_CONTENTLABEL_H
