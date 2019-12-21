/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace votca {
namespace tools {

using namespace std;
using namespace boost;

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
    auto it = int_vals.find(key);
    int_string_id.append(lexical_cast<string>(it->second));
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
    double_string_id.append(sig_fig_(it->second, 8));
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
    str_string_id.append(lexical_cast<string>(it->second));
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
  int_vals_ = int_vals;
  initStringId_();
}

void GraphNode::setDouble(const unordered_map<string, double> double_vals) {
  double_vals_ = double_vals;
  initStringId_();
}

void GraphNode::setStr(const unordered_map<string, string> str_vals) {
  str_vals_ = str_vals;
  initStringId_();
}

Index GraphNode::getInt(const string str) {
  if (int_vals_.count(str) == 0) {
    throw invalid_argument(
        "GraphNode does not "
        "contain value");
  }
  return int_vals_[str];
}

double GraphNode::getDouble(const string str) {
  if (double_vals_.count(str) == 0) {
    throw invalid_argument(
        "GraphNode does not "
        "contain value");
  }
  return double_vals_[str];
}

std::string GraphNode::getStr(const string str) {
  if (str_vals_.count(str) == 0) {
    throw invalid_argument(
        "GraphNode does not "
        "contain value");
  }
  return str_vals_[str];
}

bool GraphNode::operator!=(const GraphNode gn) const {
  return (str_id_.compare(gn.str_id_) != 0);
}

bool GraphNode::operator==(const GraphNode gn) const {
  return !((*this) != gn);
}

ostream& operator<<(ostream& os, const GraphNode gn) {
  os << "Integer Values" << endl;
  for (const auto& int_val : gn.int_vals_) {
    os << int_val.first << " " << int_val.second << endl;
  }
  os << "Double  Values" << endl;
  for (const auto& double_val : gn.double_vals_) {
    os << double_val.first << " " << double_val.second << endl;
  }
  os << "String  Values" << endl;
  for (const auto& str_val : gn.str_vals_) {
    os << str_val.first << " " << str_val.second << endl;
  }
  return os;
}

bool cmpNode(GraphNode gn1, GraphNode gn2) {
  string str1_Id = gn1.getStringId();
  return str1_Id.compare(gn2.getStringId()) < 0;
}
}  // namespace tools
}  // namespace votca
