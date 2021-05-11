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

  /**
   * \brief Resolves "link" attribute in the Property by filling in defaults.
   * Already existing tags are not overwritten
   */
  void ResolveLinks(tools::Property &prop) const;

  /**
   * \brief Loads default options stored in defaults_path_
   */
  Property LoadDefaults(const std::string &calculatorname) const;

  void RemoveAttributesFromUserOptions(Property &user_options) const;

  /**
   * \brief Updates default options with user options
   */
  void UpdateWithUserOptions(Property &options, const Property &user_options,
                             const std::string &calcname) const;

  /**
   * \brief Updates default options with user options from commandline
   */
  void UpdateWithUserCommandline(Property &options,
                                 const std::vector<std::string> &updates) const;

  /**
   * \brief Checks that all options with default="REQUIRED" are filled in
   */
  void CheckRequired(const Property &options) const;

  /**
   * \brief Removes tags which have no value and default="OPTIONAL"
   */
  void RemoveOptional(Property &options) const;

  void InjectDefaultsAsValues(Property &options) const;

  void CheckChoices(const Property &options) const;

  void CleanAttributes(Property &options,
                       const std::vector<std::string> &attributes) const;
  /**
   * \brief Load the default options and merge them with the user input
   *
   * Defaults are overwritten with user input
   */
  Property ProcessUserInput(const Property &user_input,
                            const std::vector<std::string> &cli_input,
                            const std::string &calcname) const;

  Property CalculatorOptions(const std::string &calcname) const {
    Property print = LoadDefaults(calcname);
    ResolveLinks(print);
    CleanAttributes(print, {"link", "item"});
    return print;
  }

 private:
  void OverwriteDefaultsWithUserInput(const Property &p,
                                      Property &defaults) const;
  // Copy the defaults into the value
  static void RecursivelyCheckOptions(const Property &p);
  static bool IsValidOption(const Property &prop,
                            const std::vector<std::string> &choices);
  static std::vector<std::string> GetPropertyChoices(const Property &p);

  std::string defaults_path_;
};

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_OPTIONSHANDLER_H
