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

#include "votca/tools/attributes.h"
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
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
static string sig_fig_(double val, Index sf) {
  return ([val](Index number_of_sig_figs) -> string {
    stringstream lStream;
    lStream << setprecision(int(number_of_sig_figs)) << val;
    return lStream.str();
  })(sf);
}

static void checkKey_(const std::string& key) {
  std::vector<std::string> reserved_symbols{"=", ",", ";", "{", "}"};
  for (const std::string& symbol : reserved_symbols) {
    if (key.find(symbol) != std::string::npos) {
      std::string err_msg = "Error keys to attributes cannot contain ";
      err_msg += symbol + " symbol it is reserved for internal use.";
      throw std::invalid_argument(err_msg);
    }
  }
}

static string getLabel_(const unordered_map<string, boost::any> vals) {
  vector<string> keys;
  for (auto it : vals) {
    keys.push_back(it.first);
  }
  sort(keys.begin(), keys.end());
  std::string label = "";
  for (auto key : keys) {
    auto it = vals.find(key);
    label.append(key);
    label.append("=");
    if (it->second.type() == typeid(double)) {
      double val = boost::any_cast<double>(it->second);
      label.append(sig_fig_(val, 8));
    } else if (it->second.type() == typeid(std::string)) {
      std::string val = boost::any_cast<std::string>(it->second);
      label.append(val);
    } else if (it->second.type() == typeid(Index)) {
      Index val = boost::any_cast<Index>(it->second);
      label.append(lexical_cast<std::string>(val));
    } else {
      std::string error_msg = "Unable to compile attribute label for type ";
      error_msg +=
          string(it->second.type().name()) + " currently not supported";
      throw std::runtime_error(error_msg);
    }
    label.append(",");
  }
  return label;
}

static string getLabelBrief_(const unordered_map<string, boost::any> vals) {
  vector<string> keys;
  for (auto it : vals) {
    keys.push_back(it.first);
  }
  sort(keys.begin(), keys.end());
  std::string label = "";
  for (auto key : keys) {
    auto it = vals.find(key);
    if (it->second.type() == typeid(double)) {
      double val = boost::any_cast<double>(it->second);
      label.append(sig_fig_(val, 8));
    } else if (it->second.type() == typeid(std::string)) {
      std::string val = boost::any_cast<std::string>(it->second);
      label.append(val);
    } else if (it->second.type() == typeid(Index)) {
      Index val = boost::any_cast<Index>(it->second);
      label.append(lexical_cast<std::string>(val));
    } else {
      std::string error_msg = "Unable to compile attribute label for type ";
      error_msg +=
          string(it->second.type().name()) + " currently not supported";
      throw std::runtime_error(error_msg);
    }
    label.append(",");
  }
  return label;
}

///////////////////////////////////////////////////////////
// Private Functions
///////////////////////////////////////////////////////////
/// Used to reinitialize the string id if any of the contents
/// of the graphnode change
void Attributes::buildLabels_() {
  full_label_.clear();
  full_label_.append(getLabel_(attributes_));
  if (full_label_.length() > 0) full_label_.back() = ';';

  brief_label_.clear();
  brief_label_.append(getLabelBrief_(attributes_));
  if (brief_label_.length() > 0) brief_label_.back() = ';';
}

///////////////////////////////////////////////////////////
// Public Functions
///////////////////////////////////////////////////////////
Attributes::Attributes(const unordered_map<string, boost::any> attributes) {
  attributes_ = attributes;
  buildLabels_();
}

void Attributes::set(const unordered_map<string, boost::any> attributes) {
  for (auto& pr : attributes_) {
    checkKey_(pr.first);
  }
  attributes_ = attributes;
  buildLabels_();
}

void Attributes::reset(std::string key, const boost::any value) {
  checkKey_(key);
  if (attributes_.count(key) == 0) {
    throw std::invalid_argument("Cannot reset attribute, " + key +
                                " the key is not known.");
  }
  attributes_[key] = value;
  buildLabels_();
}

void Attributes::add(std::string key, const boost::any value) {
  assert(attributes_.count(key) == 0 &&
         "Cannot add attribute Attributes instance, the key has already been "
         "used");
  checkKey_(key);
  attributes_[key] = value;
  buildLabels_();
}

bool Attributes::exists(const std::string& key) const {
  return attributes_.count(key) > 0;
}

bool Attributes::operator!=(const Attributes attr) const {
  return (full_label_.compare(attr.full_label_) != 0);
}

bool Attributes::operator==(const Attributes attr) const {
  return !((*this) != attr);
}

ostream& operator<<(ostream& os, const Attributes attr) {
  os << "Integer Values" << endl;
  for (const auto& val : attr.attributes_) {
    if (val.second.type() == typeid(int)) {
      int v = boost::any_cast<int>(val.second);
      os << val.first << " " << v << endl;
    }
  }
  os << "Double  Values" << endl;
  for (const auto& val : attr.attributes_) {
    if (val.second.type() == typeid(double)) {
      double v = boost::any_cast<double>(val.second);
      os << val.first << " " << v << endl;
    }
  }
  os << "String  Values" << endl;
  for (const auto& val : attr.attributes_) {
    if (val.second.type() == typeid(std::string)) {
      std::string v = boost::any_cast<std::string>(val.second);
      os << val.first << " " << v << endl;
    }
  }
  return os;
}

bool cmpAttributes(Attributes attr1, Attributes attr2) {
  string label = attr1.getContentLabel();
  return label.compare(attr2.getContentLabel()) < 0;
}
}  // namespace tools
}  // namespace votca
