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

#include "votca/tools/contentlabel.h"
#include "votca/tools/types.h"

#include <algorithm>
#include <functional>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

using namespace std;
namespace votca {
namespace tools {

static size_t generateHash(vector<string> labels) {
  string flat_label = "";
  for (string label : labels) {
    flat_label.append(label);
  }
  return hash<std::string>{}(flat_label);
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
}

void ContentLabel::append_(ContentLabel label) {
  full_label_.insert(full_label_.end(), label.full_label_.begin(),
                     label.full_label_.end());
  brief_label_.insert(brief_label_.end(), label.brief_label_.begin(),
                      label.brief_label_.end());
}

bool ContentLabel::containsBranch_() {
  if (full_label_.back() == "}") return true;
  return false;
}

void ContentLabel::checkKey_(const std::string& key) {
  std::vector<std::string> reserved_symbols{"=", ",", ";", "{", "}"};
  for (const std::string& symbol : reserved_symbols) {
    if (key.find(symbol) != std::string::npos) {
      std::string err_msg = "Error keys to ContentLabel cannot contain ";
      err_msg += symbol + " symbol it is reserved for internal use.";
      throw std::invalid_argument(err_msg);
    }
  }
}
void ContentLabel::ContentLabel(std::unordered_map<string, boost::any> values) {
  full_label_ = getLabel_(values);
  if (full_label_.size() > 0) full_label_.back().back() = ';';

  brief_label_ = getLabelBrief_(values);
  if (brief_label_.size() > 0) brief_label_.back().back() = ';';

  hash_ = generateLabel(full_label_);
}

std::string get(LabelType type = LabelType::verbose) {
  if (type = LabelType::verbose) return full_label_;
  return brief_label_;
}

void ContentLabel::add(GraphNode gn) {
  if (full_label_.back().back() == "}") {
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

  return;
}
void ContentLabel::add(ContentLabel label) {
  if (containsBranch()) {
  }
  return;
}

void ContentLabel::makeBranch() {
  if (containsBranch()) return;
  full_label_.push_front("{");
  full_label_.push_back("}");
  brief_label_.push_front("{");
  brief_label_.push_back("}");
}

bool ContentLabel::operator!=(const ContentLabel label) const {
  if (label.hash_ != hash_) return true;
  if (label.labels_.size() != labels_.size()) return true;
}
bool ContentLabel::operator==(const ContentLabel label) const;
bool ContentLabel::operator<(const ContentLabel label) const;
bool ContentLabel::operator<=(const ContentLabel label) const;
bool ContentLabel::operator>(const ContentLabel label) const;
bool ContentLabel::operator>=(const ContentLabel label) const;

}  // namespace tools
}  // namespace votca
