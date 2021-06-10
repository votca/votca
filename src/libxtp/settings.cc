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

// Local VOTCA includes
#include "votca/xtp/settings.h"

namespace votca {
namespace xtp {

void Settings::read_property(const votca::tools::Property& properties,
                             const std::string& key) {
  for (const auto& prop : properties.get(key)) {
    const auto& name = prop.name();
    this->nodes_[name] = prop;
  }
}

void Settings::load_from_xml(const std::string& path) {
  votca::tools::Property options;
  options.LoadFromXML(path);
  this->read_property(options, root_key_);
}

// TODO: move this method to Property
void Settings::amend(const Settings& other) {
  // Merge general properties
  for (const auto& pair : other.nodes_) {
    auto it = this->nodes_.find(pair.first);
    if (it == this->nodes_.end()) {
      const auto& val = pair.second;
      this->nodes_[pair.first] = val;
    }
  }
}

bool Settings::has_key(const std::string& key) const {
  auto it = this->nodes_.find(key);
  return (it != this->nodes_.end()) ? true : false;
}

void Settings::add(const std::string& key, const std::string& value) {
  std::string primary_key = key.substr(0, key.find("."));
  std::string secondary_key = key.substr(key.find(".") + 1);
  votca::tools::Property& prop = this->nodes_[primary_key];
  prop.add(secondary_key, value);
}

void Settings::validate() const {
  std::vector<std::string> keywords = mandatory_keyword_;
  if (this->get("name") != "xtp") {
    this->check_mandatory_keyword("executable");
  }
  for (const auto& x : keywords) {
    this->check_mandatory_keyword(x);
  }
  std::stringstream stream;
  // Check that the input keys are valid
  for (const auto& pair : this->nodes_) {
    auto it = std::find(this->general_properties_.cbegin(),
                        this->general_properties_.cend(), pair.first);

    if (it == this->general_properties_.cend()) {
      stream << "Unrecognized keyword \"" << pair.first << "\"\n"
             << "Keywords must be one of the following:\n";
      for (const std::string& key : this->general_properties_) {
        stream << key << "\n";
      }
      throw std::runtime_error(stream.str());
    }
  }
}

void Settings::check_mandatory_keyword(const std::string& key) const {
  std::stringstream stream;
  auto it = this->nodes_.find(key);
  if (it == this->nodes_.end()) {
    stream << "the " << key << " keyword is mandatory\n";
    auto it2 = this->keyword_options_.find(key);
    if (it2 != this->keyword_options_.end()) {
      stream << key << "must be one of the following values:\n";
      for (const auto& x : it2->second) {
        stream << x << "\n";
      }
    }
    throw std::runtime_error(stream.str());
  }
}

std::ostream& operator<<(std::ostream& out, const Settings& sett) {
  for (const auto& pair : sett.nodes_) {
    out << "key: " << pair.first << "\n"
        << "value: " << pair.second << "\n";
  }
  return out;
}

tools::Property Settings::to_property(const std::string& name) const {
  tools::Property root{"options", "", "."};
  tools::Property& prop = root.add(name, "");
  for (const auto& pair : this->nodes_) {
    tools::Property& new_prop = prop.add(pair.first, "");
    new_prop = pair.second;
  }
  return root;
}

}  // namespace xtp
}  // namespace votca
