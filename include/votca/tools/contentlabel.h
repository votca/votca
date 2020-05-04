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

#include <deque>
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
 private:
  std::deque<std::string> labels_;
  size_t hash_ = 0;
  size_t label_char_len_ = 0;
  LabelType type_;
  void append_(ContentLabel label);

 public:
  ContentLabel(std::unordered_map<std::string, boost::any> values,
               LabelType type);
  void add(GraphNode gn);
  void add(Branch br);
  void add(ContentLabel);

  void makeBranch();

  LabelType getLabelType() const noexcept { return type_; }

  std::string get() const;

  bool operator!=(const ContentLabel label) const;
  bool operator==(const ContentLabel label) const;
  bool operator<(const ContentLabel label) const;
  bool operator<=(const ContentLabel label) const;
  bool operator>(const ContentLabel label) const;
  bool operator>=(const ContentLabel label) const;
};

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_CONTENTLABEL_H
