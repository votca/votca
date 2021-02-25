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

// Standard includes
#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

// Third party includes
#include <boost/lexical_cast.hpp>

// Local VOTCA includes
#include "votca/tools/contentlabel.h"
#include "votca/tools/types.h"

using namespace boost;

namespace votca {
namespace tools {

  /**
   * @brief Method will generate a hash based on the contents labels
   *
   * @param labels - Each element in labels is of the KeyValType, the 
   * flat_label is constructed by appending these strings together, 
   * the actual hash is created from the flat label. 
   *
   * E.g. typically a key value pair is in the center
   *
   * key: age 
   * value: 21
   *
   * There may or may not exist a closing or opening type e.g.
   *
   * Opening: {
   * Closing: }
   *
   * The flat label then becomes
   *
   * {age=21}
   *
   * @return 
   */
static size_t generateHash(std::list<std::vector<BaseContentLabel::KeyValType>> labels) {
  std::string flat_label = "";
  for ( const auto & node : labels) {
    for (const auto & key_val_type : node) {
      flat_label = "" + key_val_type[0] + key_val_type[1] + "=" + 
        key_val_type[2] + key_val_type[3];
    }
  }
  return std::hash<std::string>{}(flat_label);
}

static bool stringsEqual(const std::string& str1, const std::string& str2) {
  return str1.compare(str2) == 0;
}

static bool stringLess(const std::string& str1, const std::string& str2) {
  if (str1.size() < str2.size()) return true;
  if (str2.size() < str1.size()) return false;
  return str1.compare(str2) > 0;
}

static void checkString_(const std::string& key) {
  std::vector<std::string> reserved_symbols{"=", ",", ";", "{", "}"};
  for (const std::string& symbol : reserved_symbols) {
    if (key.find(symbol) != std::string::npos) {
      std::string err_msg =
          "Error keys and values to BaseContentLabel cannot contain ";
      err_msg += symbol + " symbol it is reserved for internal use.";
      throw std::invalid_argument(err_msg);
    }
  }
}

static std::string sig_fig_(double val, Index sf) {
  return ([val](Index number_of_sig_figs) -> std::string {
      std::stringstream lStream;
    lStream << std::setprecision(int(number_of_sig_figs)) << val;
    return lStream.str();
  })(sf);
}

static std::string convertToString_(
    const std::unordered_map<std::string, boost::any>::iterator it) {
  if (it->second.type() == typeid(double)) {
    double val = boost::any_cast<double>(it->second);
    return sig_fig_(val, 8);
  } else if (it->second.type() == typeid(std::string)) {
    return boost::any_cast<std::string>(it->second);
  } else if (it->second.type() == typeid(Index)) {
    Index val = boost::any_cast<Index>(it->second);
    return lexical_cast<std::string>(val);
  } else if (it->second.type() == typeid(int)) {
    int val = boost::any_cast<int>(it->second);
    return lexical_cast<std::string>(val);
  }
  std::string error_msg = "Unable to compile attribute label for type ";
  error_msg += std::string(it->second.type().name()) + " currently not supported";
  throw std::runtime_error(error_msg);
}

static std::string 
buildVerboseNodeLabel_(const std::vector<BaseContentLabel::KeyValType> & node ) {
  std::string label = "";
  for (const BaseContentLabel::KeyValType& key_val_type : node) {
    label += key_val_type[0];
    label += key_val_type[1];
    label += "=";
    label += key_val_type[2];
    label += key_val_type[3];
  }
  return label;
}

static std::string buildLabel_(
    const std::list<std::vector<BaseContentLabel::KeyValType>>& values,
    const LabelType& type) {

  std::string label = "";
  if (type == LabelType::verbose) {
    for ( const auto & node : values ) {
      label += buildVerboseNodeLabel_(node);
    }
  } else if(type == LabelType::verbose_without_end_nodes){
    if( values.size() > 2) {
      auto iter_stop = values.end();
      --iter_stop;
      auto iter_start = values.begin();
      ++iter_start;
      for( auto iter = iter_start; iter != iter_stop; ++iter) {
        label += buildVerboseNodeLabel_(*iter);
      }
    }
  } else if(type == LabelType::verbose_start_node){
    label += buildVerboseNodeLabel_(*values.begin());
  } else if(type == LabelType::verbose_terminating_node){
    auto iter_stop = values.end();
    --iter_stop;
    label += buildVerboseNodeLabel_(*(iter_stop));
  } else {
    for ( const auto & node : values ) {
      for (const BaseContentLabel::KeyValType& key_val_type : node) {
        label += key_val_type[0];
        // Skip the key
        label += key_val_type[2];
        label += key_val_type[3];
      }
    }
  }
  return label;
}

/*
static vector<string> buildKeys_(const unordered_map<string, boost::any> vals) {
  vector<string> keys_temp;
  for (auto it : vals) {
    keys_temp.push_back(it.first);
  }
  sort(keys_temp.begin(), keys_temp.end());
  vector<std::string> keys;
  for (auto key_t : keys_temp) {
    auto it = vals.find(key);
    std::string key;
    key.append(key_t);
    key.append("=");
    keys.push_back(key);
  }
  return keys;
}

static vector<string> buildValues_(
    const unordered_map<string, boost::any> vals) {
  vector<string> keys;
  for (auto it : vals) {
    keys.push_back(it.first);
  }
  sort(keys.begin(), keys.end());
  vector<string> labels;
  for (auto key : keys) {
    auto it = vals.find(key);
    std::string label = "";
    label.append(convertToString_(it));
    label.append(",");
    labels.push_back(label);
  }
  return labels;
}*/

// Purpose is to place the contents in the vector
void BaseContentLabel::initLabels_(
    std::unordered_map<std::string, boost::any> values) {
  // Sort the keys alphabetically
  std::vector<std::string> keys_temp;
  for (auto it : values) {
    keys_temp.push_back(it.first);
  }
  sort(keys_temp.begin(), keys_temp.end());

  std::vector<KeyValType> node;
  // Add the labels with their corresponding keys in the correct order
  for (auto key : keys_temp) {
    auto it = values.find(key);
    std::string val = convertToString_(it);
    checkString_(key);
    checkString_(val);
    KeyValType key_val_type = {"", key, val, ","};
    node.push_back(key_val_type);
  }
  // Change the last type from a comma to a semicolon to indicate end of
  // the node
  node.back()[3] = ";";

  // Update the char length of the BaseContentLabel
  for (KeyValType& element : node) {
    label_char_len_ += element[0].size();
    label_char_len_ += element[1].size();
    label_char_len_ += element[2].size();
    label_char_len_ += element[3].size();
  }

  labels_.push_back(node);
  // Build the hash
  hash_ = generateHash(labels_);
}

BaseContentLabel::BaseContentLabel(std::unordered_map<std::string, boost::any> values) {
  initLabels_(values);
}

std::string BaseContentLabel::get(const LabelType& type) const {
  return buildLabel_(labels_, type);
}

void BaseContentLabel::clear() {
  hash_ = 0;
  label_char_len_ = 0;
  labels_.clear();
}

void BaseContentLabel::append_(std::list<std::vector<KeyValType>> label) {
  labels_.insert(labels_.end(), label.begin(), label.end());
  hash_ = generateHash(labels_);
}

void BaseContentLabel::append(BaseContentLabel label) {
  label_char_len_ += label.label_char_len_;
  append_(label.labels_);
}

void BaseContentLabel::reverse() {
  labels_.reverse();
}

bool BaseContentLabel::operator!=(const BaseContentLabel& label) const {
  if (label.hash_ != hash_) return true;
  if (label.label_char_len_ != label_char_len_) return true;
  if (label.labels_.size() != labels_.size()) return true;

  auto iter1 = labels_.begin();
  auto iter2 = label.labels_.begin();

  while (iter1 != labels_.end() ) {
    
    if( iter1->size() != iter2->size() ) {
      return false;
    } else {
      size_t num_elements_in_node = iter1->size(); 
      for (size_t ind = 0; ind < num_elements_in_node; ++ind) {
        for (size_t ind2 = 0; ind2 < arr_len; ++ind2) {
          if (not stringsEqual(iter1->at(ind)[ind2], iter2->at(ind)[ind2])) {
            return true;
          }
        }
      }
    }
    ++iter1;
    ++iter2;
  }
  return false;
}

bool BaseContentLabel::operator==(const BaseContentLabel& label) const {
  return not BaseContentLabel::operator!=(label);
}

bool BaseContentLabel::operator<(const BaseContentLabel& label) const {
  if (label_char_len_ < label.label_char_len_)
    return true;
  else if (label_char_len_ > label.label_char_len_)
    return false;
  // Which one has fewer elements
  size_t num_elements = labels_.size();
  if (label.labels_.size() < num_elements) {
    num_elements = label.labels_.size();
  }

  auto iter1 = labels_.begin();
  auto iter2 = label.labels_.begin();

  while (iter1 != labels_.end() ) {
    
    if( iter1->size() > iter2->size() ) {
      return false;
    } else if(iter1->size() < iter2->size() ) {
      return true;  
    } else {
      size_t num_elements_in_node = iter1->size(); 
      for (size_t ind = 0; ind < num_elements_in_node; ++ind) {
        for (size_t ind2 = 0; ind2 < arr_len; ++ind2) {
          if (stringLess(iter1->at(ind)[ind2], iter2->at(ind)[ind2])) {
            return true;
          }
        }
      }
    }
    ++iter1;
    ++iter2;
  }

  return false;
}

bool BaseContentLabel::operator<=(const BaseContentLabel& label) const {
  if (BaseContentLabel::operator==(label)) return true;
  return BaseContentLabel::operator<(label);
}

bool BaseContentLabel::operator>(const BaseContentLabel& label) const {
  return not BaseContentLabel::operator<=(label);
}

bool BaseContentLabel::operator>=(const BaseContentLabel& label) const {
  return not BaseContentLabel::operator<(label);
}

}  // namespace tools
}  // namespace votca
