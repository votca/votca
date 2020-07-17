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

#ifndef VOTCA_XTP_SETTINGS_H
#define VOTCA_XTP_SETTINGS_H

// Standard includes
#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>

// VOTCA includes
#include <votca/tools/property.h>

/*
 * \brief Hierarchical representation of a QM Package input
 */
namespace votca {
namespace xtp {

class Settings {
 public:
  // Decompose a votca::tools::Property object into Settings
  Settings() = default;
  Settings(const std::string& root_key) : _root_key{root_key} {};
  ~Settings() = default;

  /**
   * \brief Transform Properties into settings
   * @param Property object
   */
  void read_property(const votca::tools::Property& properties,
                     const std::string& key);

  /**
   * \brief read Settings from xml
   * @param path to the xml
   */
  void load_from_xml(const std::string& path);

  /**
   * \brief Fill the missing values using another Settings instance
   * @param other Settings object
   */
  void amend(const Settings& other);

  /**
   * \brief Add property
   * @param key to the property
   * @param value content of the property
   */
  void add(const std::string& key, const std::string& value);

  /**
   * \brief Get a given key
   * @param key
   */
  template <typename T = std::string>
  T get(const std::string& key) const {
    auto delimiter = ".";
    std::string primary_key = key;
    std::string secondary_key;

    auto iter = key.find(delimiter);
    if (iter != std::string::npos) {
      primary_key = key.substr(0, key.find(delimiter));
      secondary_key = key.substr(key.find(delimiter) + 1);
    }

    auto it = this->_nodes.find(primary_key);
    if (it == this->_nodes.end()) {
      std::ostringstream oss;
      oss << "Unknown keyword: " << key << "\n";
      throw std::runtime_error(oss.str());
    } else if (secondary_key.empty()) {
      return it->second.as<T>();
    } else {
      return it->second.get(secondary_key).as<T>();
    }
  }
  /**
   * \brief Retrieve property
   * @param key to property
   */
  const votca::tools::Property& property(const std::string& key) const {
    return _nodes.at(key);
  }

  /**
   * \brief Check if a property exists
   * @param key
   */
  bool has_key(const std::string& key) const;

  /**
   * \brief Check that the input is correct
   */
  void validate() const;

  /**
   * \brief Convert a Settings object into a Property
   * @param root name
   */
  tools::Property to_property(const std::string& name) const;

  friend std::ostream& operator<<(std::ostream& out, const Settings& sett);

 private:
  using Settings_map = std::unordered_map<std::string, votca::tools::Property>;
  Settings_map _nodes;  // {Key, Property} Map
  std::string _root_key;

  std::string get_primary_key(const std::string& key) {
    return key.substr(0, key.find("."));
  }

  void check_mandatory_keyword(const std::string& key) const;

  std::vector<std::string> _general_properties = {
      "auxbasisset",            // string
      "basisset",               // string
      "charge",                 // index
      "cleanup",                // string
      "convergence_tightness",  // std::string
      "dipole_spacing",         // boolean
      "ecp",                    // string
      "executable",             // string
      // "external_charge",        // Eigen::Vector9d
      "functional",     // string
      "name",           // string
      "optimize",       // boolean
      "orca",           // string
      "polarization",   // boolean
      "read_guess",     // boolean
      "spin",           // index
      "scratch",        // string
      "write_charges",  // boolean
      "xtpdft",         // string
  };

  std::vector<std::string> _mandatory_keyword = {
      "functional",  // string
      "name",        // string, one of: orca
  };

  std::unordered_map<std::string, std::vector<std::string>> _keyword_options{
      {"name", {"orca"}}};
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_SETTINGS_H
