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

#include "../../include/votca/tools/graphnode.h"
#include <iostream>
/*#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
*/
namespace votca {
namespace tools {

using namespace std;
/*using namespace boost;
///////////////////////////////////////////////////////////
// Local Functions
///////////////////////////////////////////////////////////
/// Converts a double into a string with max number of significant
/// figures indicated by sf
string sig_fig_(double val, Index sf) {
  return ([val](Index number_of_sig_figs) -> string {
    stringstream lStream;
    lStream << setprecision(int(number_of_sig_figs)) << val;
    return lStream.str();
  })(sf);
}

void checkLabel_(const std::string& label) {
  if (label.find("=") != std::string::npos)
    throw std::invalid_argument(
        "Error labels to graphnodes cannot contain = symbol it is reserved for "
        "internal use.");
  if (label.find(",") != std::string::npos)
    throw std::invalid_argument(
        "Error labels to graphnodes cannot contain , symbol it is reserved for "
        "internal use.");
  if (label.find(";") != std::string::npos)
    throw std::invalid_argument(
        "Error labels to graphnodes cannot contain ; symbol it is reserved for "
        "internal use.");
}

/// Turns the map of ints into a string that is sorted alphabetically
/// by the keys
string getIntStringId_(const unordered_map<string, Index> int_vals) {
  vector<string> keys;
  // Grab integer keys sort alphabetically and store in string_id
  string int_string_id;
  for (auto it : int_vals) {
    keys.push_back(it.first);
  }
  sort(keys.begin(), keys.end());
  for (auto key : keys) {
    int_string_id.append(key);
    int_string_id.append("=");
    auto it = int_vals.find(key);
    int_string_id.append(lexical_cast<string>(it->second));
    int_string_id.append(",");
  }
  return int_string_id;
}

/// Turns the map of ints into a string that is sorted alphabetically
/// by the keys, but ignores the keys in actually building the id
string getIntStringIdIgnoreLabel_(const unordered_map<string, Index> int_vals) {
  vector<string> keys;
  // Grab integer keys sort alphabetically and store in string_id
  string int_string_id;
  for (auto it : int_vals) {
    keys.push_back(it.first);
  }
  sort(keys.begin(), keys.end());
  for (auto key : keys) {
    auto it = int_vals.find(key);
    int_string_id.append(lexical_cast<string>(it->second));
    int_string_id.append(",");
  }
  return int_string_id;
}

/// Turns the map of doubles into a string that is sorted alphabetically
/// by the keys
string getDoubleStringId_(const unordered_map<string, double> double_vals) {
  vector<string> keys;
  // Grab double keys sort alphabetically and store in string_id
  string double_string_id;
  for (auto it : double_vals) {
    keys.push_back(it.first);
  }
  sort(keys.begin(), keys.end());
  for (auto key : keys) {
    double_string_id.append(key);
    auto it = double_vals.find(key);
    double_string_id.append("=");
    double_string_id.append(sig_fig_(it->second, 8));
    double_string_id.append(",");
  }
  return double_string_id;
}

/// Turns the map of doubles into a string that is sorted alphabetically
/// by the keys, but ignore label
string getDoubleStringIdIgnoreLabel_(
    const unordered_map<string, double> double_vals) {
  vector<string> keys;
  // Grab double keys sort alphabetically and store in string_id
  string double_string_id;
  for (auto it : double_vals) {
    keys.push_back(it.first);
  }
  sort(keys.begin(), keys.end());
  for (auto key : keys) {
    auto it = double_vals.find(key);
    double_string_id.append(sig_fig_(it->second, 8));
    double_string_id.append(",");
  }
  return double_string_id;
}

/// Turns the map of strings into a string that is sorted alphabetically
/// by the keys
string getStrStringId_(const unordered_map<string, string> str_vals) {
  vector<string> keys;
  // Grab string keys sort alphabetically and store in string_id
  string str_string_id;
  for (auto it : str_vals) {
    keys.push_back(it.first);
  }
  sort(keys.begin(), keys.end());
  for (auto key : keys) {
    str_string_id.append(key);
    auto it = str_vals.find(key);
    str_string_id.append("=");
    str_string_id.append(lexical_cast<string>(it->second));
    str_string_id.append(",");
  }
  return str_string_id;
}

/// Turns the map of strings into a string that is sorted alphabetically
/// by the keys, but ignores the key in building the id
string getStrStringIdIgnoreLabel_(
    const unordered_map<string, string> str_vals) {
  vector<string> keys;
  // Grab string keys sort alphabetically and store in string_id
  string str_string_id;
  for (auto it : str_vals) {
    keys.push_back(it.first);
  }
  sort(keys.begin(), keys.end());
  for (auto key : keys) {
    auto it = str_vals.find(key);
    str_string_id.append(lexical_cast<string>(it->second));
    str_string_id.append(",");
  }
  return str_string_id;
}

///////////////////////////////////////////////////////////
// Private Functions
///////////////////////////////////////////////////////////
/// Used to reinitialize the string id if any of the contents
/// of the graphnode change
void GraphNode::initStringId_() {
  str_id_.clear();
  str_id_.append(getIntStringId_(int_vals_));
  str_id_.append(getDoubleStringId_(double_vals_));
  str_id_.append(getStrStringId_(str_vals_));
  if (str_id_.length() > 0) str_id_.back() = ';';
  str_id_no_label_.clear();
  str_id_no_label_.append(getIntStringIdIgnoreLabel_(int_vals_));
  str_id_no_label_.append(getDoubleStringIdIgnoreLabel_(double_vals_));
  str_id_no_label_.append(getStrStringIdIgnoreLabel_(str_vals_));
  if (str_id_no_label_.length() > 0) str_id_no_label_.back() = ';';
}

///////////////////////////////////////////////////////////
// Public Functions
///////////////////////////////////////////////////////////
GraphNode::GraphNode(const unordered_map<string, Index> int_vals,
                     const unordered_map<string, double> double_vals,
                     const unordered_map<string, string> str_vals) {
  int_vals_ = int_vals;
  double_vals_ = double_vals;
  str_vals_ = str_vals;
  initStringId_();
}

void GraphNode::setInt(const unordered_map<string, Index> int_vals) {
  for (auto& pr : int_vals) {
    checkLabel_(pr.first);
  }
  int_vals_ = int_vals;
  initStringId_();
}

void GraphNode::setDouble(const unordered_map<string, double> double_vals) {
  for (auto& pr : double_vals) {
    checkLabel_(pr.first);
  }
  double_vals_ = double_vals;
  initStringId_();
}

void GraphNode::setStr(const unordered_map<string, string> str_vals) {
  for (auto& pr : str_vals) {
    checkLabel_(pr.first);
  }
  str_vals_ = str_vals;
  initStringId_();
}

void GraphNode::resetInt(std::string label, const Index& value) {
  checkLabel_(label);
  if (int_vals_.count(label) == 0) {
    throw std::invalid_argument(
        "Cannot reset graph node int, the label is not known.");
  }
  int_vals_[label] = value;
  initStringId_();
}

void GraphNode::resetDouble(std::string label, const double& value) {
  checkLabel_(label);
  if (double_vals_.count(label) == 0) {
    throw std::invalid_argument(
        "Cannot reset graph node double, the label is not known.");
  }
  double_vals_[label] = value;
  initStringId_();
}

void GraphNode::resetStr(std::string label, const std::string& value) {
  checkLabel_(label);
  if (str_vals_.count(label) == 0) {
    throw std::invalid_argument(
        "Cannot reset graph node string, the label is not known.");
  }
  str_vals_[label] = value;
  initStringId_();
}

void GraphNode::addInt(std::string label, const Index value) {
  assert(int_vals_.count(label) == 0 &&
         "Cannot add int to GraphNode label has already been used");
  checkLabel_(label);
  int_vals_[label] = value;
  initStringId_();
}

void GraphNode::addDouble(std::string label, const double value) {
  assert(double_vals_.count(label) == 0 &&
         "Cannot add double to GraphNode label has already been used");
  checkLabel_(label);
  double_vals_[label] = value;
  initStringId_();
}

void GraphNode::addStr(std::string label, const std::string& value) {
  assert(str_vals_.count(label) == 0 &&
         "Cannot add string to GraphNode label has already been used");
  checkLabel_(label);
  str_vals_[label] = value;
  initStringId_();
}

bool GraphNode::hasInt(const std::string& label) const {
  return int_vals_.count(label) > 0;
}

Index GraphNode::getInt(const string& str) const {
  if (int_vals_.count(str) == 0) {
    throw invalid_argument(
        "GraphNode does not "
        "contain value");
  }
  return int_vals_.at(str);
}

double GraphNode::getDouble(const string& str) const {
  if (double_vals_.count(str) == 0) {
    throw invalid_argument(
        "GraphNode does not "
        "contain value");
  }
  return double_vals_.at(str);
}

std::string GraphNode::getStr(const string& str) const {
  if (str_vals_.count(str) == 0) {
    throw invalid_argument(
        "GraphNode does not "
        "contain value");
  }
  return str_vals_.at(str);
}

bool GraphNode::operator!=(const GraphNode gn) const {
  return (str_id_.compare(gn.str_id_) != 0);
}

bool GraphNode::operator==(const GraphNode gn) const {
  return !((*this) != gn);
}
*/
ostream& operator<<(ostream& os, const GraphNode gn) {
  os << "Integer Values" << endl;
  std::unordered_map<std::string, int> int_vals = gn.getAll<int>();
  for (const auto& int_val : int_vals) {
    os << int_val.first << " " << int_val.second << endl;
  }
  os << "Double  Values" << endl;
  std::unordered_map<std::string, double> double_vals = gn.getAll<double>();
  for (const auto& double_val : double_vals) {
    os << double_val.first << " " << double_val.second << endl;
  }
  os << "String  Values" << endl;
  std::unordered_map<std::string, std::string> str_vals =
      gn.getAll<std::string>();
  for (const auto& str_val : str_vals) {
    os << str_val.first << " " << str_val.second << endl;
  }
  return os;
}

bool cmpNode(GraphNode gn1, GraphNode gn2) {
  ContentLabel label = gn1.getContentLabel();
  return gn2.getContentLabel() < label;
}
}  // namespace tools
}  // namespace votca*/
