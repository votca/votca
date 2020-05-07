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

#ifndef VOTCA_TOOLS_CALCULATOR_H
#define VOTCA_TOOLS_CALCULATOR_H

#include "globals.h"
#include "property.h"
#include "propertyiomanipulator.h"

namespace votca {
namespace tools {

/**
 * \brief Base class for all calculators
 *
 * Calculators are grouped in CalculatorFactories and are run by Threads
 * or Applications. Every calculator has a description (an XML file) installed
 * in VOTCASHARE which is used to compile HELP.
 * This XML file also contains default values.
 *
 */
class Calculator {
 public:
  Calculator() = default;
  virtual ~Calculator() = default;
  /**
   * \brief Calculator name
   *
   * This name is used to register a calculator in a Factory
   * It the name of the XML file with the default calculator options
   * stored in VOTCASHARE
   *
   * @return calculator name
   */
  virtual std::string Identify() = 0;
  /**
   * \brief Initializes a calculator from an XML file with options
   *
   * Options are passed to a calculator by the Application
   * These option overwrite defaults
   *
   * @param options Property object passed by the application to a calculator
   */
  virtual void Initialize(const Property &user_options) = 0;
  /**
   * \brief Sets number of threads to use
   *
   * If only one thread is used, this calculator behaves as a master
   *
   * @param nThreads number of threads running this calculator
   *
   */
  void setnThreads(Index nThreads) {
    _nThreads = nThreads;
    _maverick = (_nThreads == 1) ? true : false;
  }
  /**
   * \brief Outputs all options of a calculator
   *
   * @param output stream
   */
  void DisplayOptions(std::ostream &out);

  /**
   * \brief Loads default options stored in VOTCASHARE
   */
  Property LoadDefaults(const std::string package = "tools");

  /**
   * \brief Updates user options with default options stored in VOTCASHARE
   *
   * If a value is not given or tag is not present and at the same time
   * a default value exists in the corresponding XML file in VOTCASHARE
   * a tag is created and/or a default value is assigned to it
   */
  void UpdateWithUserOptions(Property &default_options,
                             const Property &user_options);

  /**
   * \brief Load the default options and merge them with the user input
   *
   * Defaults are overwritten with user input
   */
  Property LoadDefaultsAndUpdateWithUserOptions(const std::string package,
                                                const Property &user_options) {
    Property defaults = LoadDefaults(package);
    InjectDefaultsAsValues(defaults);
    UpdateWithUserOptions(defaults, user_options);
    RecursivelyCheckOptions(defaults);
    return defaults;
  }

 protected:
  Index _nThreads;
  bool _maverick;

  void OverwriteDefaultsWithUserInput(const Property &p, Property &defaults);
  // Copy the defaults into the value
  static void InjectDefaultsAsValues(Property &defaults);
  static void RecursivelyCheckOptions(const Property &prop);
  static bool IsValidOption(const Property &p,
                            const std::vector<std::string> &choices);
  static std::string GetVotcaShare();
  static std::vector<std::string> GetPropertyChoices(const Property &prop);

  template <typename T>
  static bool IsValidCast(const tools::Property &prop) {
    try {
      prop.as<T>();
      return true;
    } catch (const std::runtime_error &e) {
      return false;
    }
  }
};

}  // namespace tools
}  // namespace votca

#endif /* VOTCA_TOOLS_CALCULATOR_H */
