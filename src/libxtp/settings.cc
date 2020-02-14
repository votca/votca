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

#include <votca/xtp/settings.h>

using votca::tools::Property;

namespace votca {
namespace xtp {

void Settings::read_property(const Property& properties) {
  for (const auto& prop : properties) {
    const auto& name = prop.name();
    this->_nodes[name] = prop;
  }
}

void Settings::merge(const Settings& other) {
  // Merge general properties
  for (const auto& pair : other._nodes) {
    auto it = this->_nodes.find(pair.first);
    if (it == this->_nodes.end()) {
      const auto& val = pair.second;
      this->_nodes[pair.first] = val;
    }
  }
}

Settings::Settings_map::const_iterator Settings::search_for_mandatory_keyword(
    std::string key) const {
  std::ostringstream oss;
  auto it = this->_nodes.find(key);
  oss << "the " << key << " keyword is mandatory\n";
  if (it == this->_nodes.end()) {
    throw std::runtime_error(oss.str());
  } else {
    return it;
  }
}

void Settings::validate() const {
  this->validate_name();
  this->search_for_mandatory_keyword("executable");

  std::ostringstream oss;
  for (const auto& pair : this->_nodes) {
    auto it = find(_general_properties.cbegin(), _general_properties.cend(),
                   pair.first);
    if (it == _general_properties.cend()) {
      oss << "Unknown keyword: " << pair.first << "\n";
      throw std::runtime_error(oss.str());
    }
  }
}

void Settings::validate_name() const {
  auto it = this->search_for_mandatory_keyword("name");
  std::string name = it->second.value();
  std::vector<std::string> names = {"orca", "gaussian", "nwchem", "xtpdft"};
  auto it2 = find(names.cbegin(), names.cend(), name);
  if (it2 == names.cend()) {
    throw std::runtime_error(
        "name must be one of: orca, gaussian, nwchem or xtpdft\n");
  }
}

std::ostream& operator<<(std::ostream& os, const Settings& sett) {
  for (const auto& prop : sett._nodes) {
    os << "name: " << prop.first << " value: " << prop.second.value() << "\n";
  }
  return os;
}

}  // namespace xtp
}  // namespace votca
