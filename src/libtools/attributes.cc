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

#include "../../include/votca/tools/attributes.h"
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

using namespace std;
using namespace boost;

namespace votca {
namespace tools {

std::vector<std::string> Attributes::getKeys() const noexcept {
  std::vector<std::string> keys;
  for (const auto& key_val : attributes_) {
    keys.push_back(key_val.first);
  }
  return keys;
}

std::unordered_map<std::string, std::string> Attributes::getTypes() const
    noexcept {
  std::unordered_map<std::string, std::string> keys_types;
  for (const auto& key_val : attributes_) {
    keys_types[key_val.first] = key_val.second.type().name();
  }
  return keys_types;
}

void Attributes::remove(const std::string& key) {
  assert(exists(key) && "Cannot remove item with provided key does not exist");
  attributes_.erase(key);
}

bool Attributes::exists(const std::string& key) const {
  return attributes_.count(key) > 0;
}

bool Attributes::operator!=(const Attributes attr) const {
  return label_ != attr.label_;
}

bool Attributes::operator==(const Attributes attr) const {
  return not((*this) != attr);
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
  ContentLabel label = attr1.getContentLabel();
  return label < attr2.getContentLabel();
}

}  // namespace tools
}  // namespace votca
