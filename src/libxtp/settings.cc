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

void Settings::read_property(const Property& properties,
                             const std::string& key) {
  for (const auto& prop : properties.get(key)) {
    const auto& name = prop.name();
    this->_nodes[name] = prop;
  }
}

void Settings::load_from_xml(const std::string& path) {
  Property options;
  options.LoadFromXML(path);
  this->read_property(options, _root_key);
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

bool Settings::exists(const std::string& key) const {
  auto it = this->_nodes.find(key);
  return (it != this->_nodes.end()) ? true : false;
}

Settings::Settings_map::const_iterator Settings::search_for_mandatory_keyword(
    const std::string& key) const {
  std::ostringstream oss;
  auto it = this->_nodes.find(key);
  if (it == this->_nodes.end()) {
    oss << "the " << key << " keyword is mandatory\n";
    auto it2 = this->_keyword_options.find(key);
    if (it2 != this->_keyword_options.end()) {
      oss << key << "must be one of the following values:\n";
      for (const auto& x : it2->second) {
        oss << x << "\n";
      }
    }
    throw std::runtime_error(oss.str());
  } else {
    return it;
  }
}

void Settings::validate() const {
  for (const auto& x : _mandatory_keyword) {
    this->search_for_mandatory_keyword(x);
  }

  /* std::ostringstream oss;
  for (const auto& pair : this->_nodes) {
    auto it = find(_general_properties.cbegin(), _general_properties.cend(),
                   pair.first);
    if (it == _general_properties.cend()) {
      oss << "Unknown keyword: " << pair.first << "\n";
      throw std::runtime_error(oss.str());
    }
   }*/
}

std::string Settings::create_orca_section(const std::string& key) const {
  std::ostringstream oss;
  std::string section = key.substr(key.find(".") + 1);
  oss << "%" << section << " " << this->get(key) << "\n"
      << "    end";
  return oss.str();
}


std::ostream& operator<<(std::ostream& os, const Settings& sett) {
  for (const auto& prop : sett._nodes) {
    os << "name: " << prop.first << " value: " << prop.second.value() << "\n";
  }
  return os;
}

}  // namespace xtp
}  // namespace votca
