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

namespace votca {
namespace xtp {

void Settings::read_property(const votca::tools::Property& properties,
                             const std::string& key) {
  for (const auto& prop : properties.get(key)) {
    const auto& name = prop.name();
    this->_nodes[name] = prop;
  }
}

void Settings::load_from_xml(const std::string& path) {
  votca::tools::Property options;
  options.LoadFromXML(path);
  this->read_property(options, _root_key);
}

// TODO: move this method to Property
void Settings::amend(const Settings& other) {
  // Merge general properties
  for (const auto& pair : other._nodes) {
    auto it = this->_nodes.find(pair.first);
    if (it == this->_nodes.end()) {
      const auto& val = pair.second;
      this->_nodes[pair.first] = val;
    }
  }
}

bool Settings::has_key(const std::string& key) const {
  auto it = this->_nodes.find(key);
  return (it != this->_nodes.end()) ? true : false;
}

void Settings::add(const std::string& key, const std::string& value) {
  std::string primary_key = key.substr(0, key.find("."));
  votca::tools::Property& prop = this->_nodes[primary_key];
  prop.add(key, value);
}

void Settings::validate() const {
  std::vector<std::string> keywords = _mandatory_keyword;
  if (this->get("name") != "xtp") {
    this->check_mandatory_keyword("executable");
  }
  for (const auto& x : keywords) {
    this->check_mandatory_keyword(x);
  }
  std::stringstream stream;
  // Check that the input keys are valid
  for (const auto& pair : this->_nodes) {
    auto it = std::find(this->_general_properties.cbegin(),
                        this->_general_properties.cend(), pair.first);

    if (it == this->_general_properties.cend()) {
      stream << "Unrecognized keyword \"" << pair.first << "\"\n"
             << "Keywords must be one of the following:\n";
      for (const std::string& key : this->_general_properties) {
        stream << key << "\n";
      }
      throw std::runtime_error(stream.str());
    }
  }
}

void Settings::check_mandatory_keyword(const std::string& key) const {
  std::stringstream stream;
  auto it = this->_nodes.find(key);
  if (it == this->_nodes.end()) {
    stream << "the " << key << " keyword is mandatory\n";
    auto it2 = this->_keyword_options.find(key);
    if (it2 != this->_keyword_options.end()) {
      stream << key << "must be one of the following values:\n";
      for (const auto& x : it2->second) {
        stream << x << "\n";
      }
    }
    throw std::runtime_error(stream.str());
  }
}

}  // namespace xtp
}  // namespace votca
