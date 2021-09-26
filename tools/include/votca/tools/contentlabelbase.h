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

#ifndef VOTCA_TOOLS_BASECONTENTLABEL_H
#define VOTCA_TOOLS_BASECONTENTLABEL_H
#pragma once

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/any.hpp>
namespace votca {
namespace tools {

enum class LabelType { 
  verbose,
  verbose_without_end_nodes,
  verbose_start_node,
  verbose_terminating_node,
  concise 
};

/**
 * @brief Content Label is meant to be a unique identifier that contains
 * the contents of some object
 *
 * The content label encodes all of the structural information about an
 * object into a single string.  The idea behind the content label is to
 * generate a label that is unique to the structure of a given graph object,
 * such that the objects structure could be simply rebuilt with the content
 * label alone.  Store the contents as a vector of strings that appear in
 * the correct sequence
 */
class BaseContentLabel {
 public:
  static constexpr int arr_len = 4;
  // index 0 - opening type, 1 key, 2 - val, 3 - closing type
  typedef std::array<std::string, arr_len> KeyValType;

 protected:

  /// Every time an existing content label is appended it appears as a
  /// new node in the list 
  std::list<std::vector<KeyValType>> labels_;
  size_t hash_ = 0;

  /// This member is in charge of recording the actual size of the string when
  /// the label is built
  size_t label_char_len_ = 0;

  void initLabels_(std::unordered_map<std::string, boost::any> values);

  void append_(std::list<std::vector<KeyValType>> label);


 public:
  BaseContentLabel() = default;
  BaseContentLabel(std::unordered_map<std::string, boost::any> values);
  virtual void append(BaseContentLabel);



  size_t getCharLen() const noexcept { return label_char_len_; }

  bool isEmpty() const noexcept { return labels_.size() == 0; }
  void clear();

  // Reverse the order of the sequence
  virtual void reverse();

  std::string get(const LabelType& label_type = LabelType::verbose) const;

  bool operator!=(const BaseContentLabel& label) const;
  bool operator==(const BaseContentLabel& label) const;
  bool operator<(const BaseContentLabel& label) const;
  bool operator<=(const BaseContentLabel& label) const;
  bool operator>(const BaseContentLabel& label) const;
  bool operator>=(const BaseContentLabel& label) const;
};

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_BASECONTENTLABEL_H
