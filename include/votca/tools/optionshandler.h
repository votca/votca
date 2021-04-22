/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_TOOLS_OPTIONSHANDLER_H
#define VOTCA_TOOLS_OPTIONSHANDLER_H

// Local VOTCA includes
#include "property.h"

namespace votca {
namespace tools {

class OptionsHandler {
 public:
  OptionsHandler(std::string defaults_path) : defaults_path_(defaults_path) {}

  void ResolveLinks(tools::Property &prop) const;

  /**
   * \brief Loads default options stored in defaults_path_
   */
  Property LoadDefaults(const std::string &calculatorname) const;

  /**
   * \brief Updates user options with default options stored in defaults_path_
   *
   * If a value is not given or tag is not present and at the same time
   * a default value exists in the corresponding XML file in defaults_path_
   * a tag is created and/or a default value is assigned to it
   */
  void UpdateWithUserOptions(Property &default_options,
                             const Property &user_options, const std::string& calcname) const;

  /**
   * \brief Load the default options and merge them with the user input
   *
   * Defaults are overwritten with user input
   */
  Property LoadDefaultsAndUpdateWithUserOptions(
      const Property &user_options, const std::string &calcname) const;

 private:
  void OverwriteDefaultsWithUserInput(const Property &p,
                                      Property &defaults) const;
  // Copy the defaults into the value
  static void InjectDefaultsAsValues(Property &defaults);
  static void RecursivelyCheckOptions(const Property &p);
  static bool IsValidOption(const Property &prop,
                            const std::vector<std::string> &choices);
  static std::vector<std::string> GetPropertyChoices(const Property &p);

  template <typename T>
  static bool IsValidCast(const tools::Property &prop) {
    try {
      prop.as<T>();
      return true;
    } catch (const std::runtime_error &e) {
      return false;
    }
  }

  std::string defaults_path_;
};

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_OPTIONSHANDLER_H
