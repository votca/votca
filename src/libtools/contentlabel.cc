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

#include <iostream>

#include "votca/tools/contentlabel.h"
#include "votca/tools/types.h"

#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <functional>
#include <iomanip>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

using namespace std;
using namespace boost;

namespace votca {
namespace tools {

static size_t generateHash(vector<ContentLabel::KeyValType> labels) {
  string flat_label = "";
  for (auto key_val_type : labels) {
    for (auto val : key_val_type) {
      // for ( size_t ind = 0; ind < ContentLabel::arr_len; ++ind ){
      flat_label.append(val);
    }
  }
  return hash<std::string>{}(flat_label);
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
          "Error keys and values to ContentLabel cannot contain ";
      err_msg += symbol + " symbol it is reserved for internal use.";
      throw std::invalid_argument(err_msg);
    }
  }
}

static string sig_fig_(double val, Index sf) {
  return ([val](Index number_of_sig_figs) -> string {
    stringstream lStream;
    lStream << setprecision(int(number_of_sig_figs)) << val;
    return lStream.str();
  })(sf);
}

static string convertToString_(
    const unordered_map<string, boost::any>::iterator it) {
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
  error_msg += string(it->second.type().name()) + " currently not supported";
  throw std::runtime_error(error_msg);
}

static std::string buildLabel_(
    const std::vector<ContentLabel::KeyValType>& values,
    const LabelType& type) {

  std::string label = "";
  if (type == LabelType::verbose) {
    for (const ContentLabel::KeyValType& key_val_type : values) {
      label += key_val_type[0];
      label += key_val_type[1];
      label += "=";
      label += key_val_type[2];
      label += key_val_type[3];
    }
  } else {
    for (const ContentLabel::KeyValType& key_val_type : values) {
      label += key_val_type[0];
      // Skip the key
      label += key_val_type[2];
      label += key_val_type[3];
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
    string key;
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
void ContentLabel::initLabels_(
    std::unordered_map<std::string, boost::any> values) {
  // Sort the keys alphabetically
  vector<string> keys_temp;
  for (auto it : values) {
    keys_temp.push_back(it.first);
  }
  sort(keys_temp.begin(), keys_temp.end());

  // Add the labels with their corresponding keys in the correct order
  for (auto key : keys_temp) {
    auto it = values.find(key);
    std::string val = convertToString_(it);
    checkString_(key);
    checkString_(val);
    KeyValType key_val_type = {"", key, val, ","};
    labels_.push_back(key_val_type);
  }
  // Change the last type from a comma to a semicolon to indicate end of
  // the node
  labels_.back()[3] = ";";

  // Update the char length of the ContentLabel
  for (KeyValType& element : labels_) {
    label_char_len_ += element[0].size();
    label_char_len_ += element[1].size();
    label_char_len_ += element[2].size();
    label_char_len_ += element[3].size();
  }

  // Build the hash
  hash_ = generateHash(labels_);
}

bool ContentLabel::containsBranch_() const {
  if (labels_.size() == 0) return false;
  if (labels_.front()[0] == "{") return true;
  return false;
}

ContentLabel::ContentLabel(std::unordered_map<string, boost::any> values) {
  initLabels_(values);
}

std::string ContentLabel::get(const LabelType& type) const {
  return buildLabel_(labels_, type);
}

void ContentLabel::clear() {
  hash_ = 0;
  label_char_len_ = 0;
  labels_.clear();
}
/*
void ContentLabel::add(GraphNode gn) {
  if (labels_.back()[3] == "}") {
    throw std::runtime_argument(
        "Can only add graph node to content label before it has been made into "
        "a branch")
  }
  ContentLabel label = gn.getContentLabel();

  full_label_.insert(full_label_.end(), label.full_label_.begin(),
                     label.full_label_.end());
  brief_label_.insert(brief_label_.end(), label.brief_label_.begin(),
                      label.brief_label_.end());
}
void ContentLabel::add(Branch br) {

  if (full_label_.back().back() != "}") {
    throw std::runtime_argument(
        "Can only add branches to content label after it has been made into a "
        "branch")
  }
  full_label_.insert(full_label_.end(), label.full_label_.begin(),
                     label.full_label_.end());
  brief_label_.insert(brief_label_.end(), label.brief_label_.begin(),
                      label.brief_label_.end());

}*/
void ContentLabel::append(ContentLabel label) {
  label_char_len_ += label.label_char_len_;
  labels_.insert(labels_.end(), label.labels_.begin(), label.labels_.end());
  hash_ = generateHash(labels_);
}

void ContentLabel::makeBranch() {
  std::cout << "calling contains branch" << std::endl;
  if (containsBranch_()) return;
  std::cout << "calling labels size" << std::endl;
  if (labels_.size() == 0) return;
  std::cout << "Assigning branch brackets" << std::endl;
  labels_.at(0)[0] = "{";
  labels_.back()[3] = "}";
  std::cout << "Done" << std::endl;
}

bool ContentLabel::operator!=(const ContentLabel& label) const {
  if (label.hash_ != hash_) return true;
  if (label.label_char_len_ != label_char_len_) return true;
  if (label.labels_.size() != labels_.size()) return true;
  for (size_t ind = 0; ind < labels_.size(); ++ind) {
    for (size_t ind2 = 0; ind2 < arr_len; ++ind2) {
      if (not stringsEqual(label.labels_.at(ind)[ind2],
                           labels_.at(ind)[ind2])) {
        return true;
      }
    }
  }
  return false;
}

bool ContentLabel::operator==(const ContentLabel& label) const {
  return not ContentLabel::operator!=(label);
}

bool ContentLabel::operator<(const ContentLabel& label) const {
  if (label_char_len_ < label.label_char_len_)
    return true;
  else if (label_char_len_ > label.label_char_len_)
    return false;
  // Which one has fewer elements
  size_t num_elements = labels_.size();
  if (label.labels_.size() < num_elements) {
    num_elements = label.labels_.size();
  }
  for (size_t ind = 0; ind < num_elements; ++ind) {
    for (size_t ind2 = 0; ind2 < arr_len; ++ind2) {
      if (stringLess(label.labels_.at(ind)[ind2], labels_.at(ind)[ind2])) {
        return true;
      }
    }
  }
  return false;
}

bool ContentLabel::operator<=(const ContentLabel& label) const {
  if (ContentLabel::operator==(label)) return true;
  return ContentLabel::operator<(label);
}

bool ContentLabel::operator>(const ContentLabel& label) const {
  return not ContentLabel::operator<=(label);
}

bool ContentLabel::operator>=(const ContentLabel& label) const {
  return not ContentLabel::operator<(label);
}

}  // namespace tools
}  // namespace votca
