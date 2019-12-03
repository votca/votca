/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <votca/tools/globals.h>
#include <votca/tools/property.h>
#include <votca/tools/propertyiomanipulator.h>

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
  virtual void Initialize(Property &options) = 0;
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
   * \brief Updates options with default options stored in VOTCASHARE
   *
   * If a value is not given or tag is not present and at the same time
   * a default value exists in the corresponding XML file in VOTCASHARE
   * a tag is created and/or a default value is assigned to it
   */
  void UpdateWithDefaults(Property &options, std::string package = "tools");

 protected:
  Index _nThreads;
  bool _maverick;

  void AddDefaults(Property &p, const Property &defaults);
};

}  // namespace tools
}  // namespace votca

#endif /* VOTCA_TOOLS_CALCULATOR_H */
