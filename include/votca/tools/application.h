/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#pragma once
#ifndef __VOTCA_APPLICATION_H
#define __VOTCA_APPLICATION_H

#include "property.h"
#include <boost/program_options.hpp>

namespace votca {
namespace tools {

class Application {
 public:
  Application();
  virtual ~Application();

  /**
   * \brief executes the program
   * \param argc argc from main
   * \param argv argv from main
   * \return return code
   */
  int Exec(int argc, char **argv);

  /**
   * \brief program name
   * \return string with program name
   *
   * overload this function to set the program name
   */
  virtual std::string ProgramName() = 0;

  /**
   * \brief version string of application
   * \return version string
   */
  virtual std::string VersionString() { return ""; }

  /**
   * \brief help text of application without version information
   * \param out ostream for output
   */
  virtual void HelpText(std::ostream &out) = 0;

  /**
   * \brief Initialize application data
   *
   * Initialize is called by run before parsing the command line.
   * All necessary command line arguments can be added here
   */
  virtual void Initialize() = 0;

  /**
   * \brief Process command line options
   * \return true to continue, false to stop
   *
   * EvaluateOptions is called by Run after parsing the command line.
   * return true if everything is ok, false to stop and show help text.
   */
  virtual bool EvaluateOptions() = 0;

  /**
   * \brief Check weather required option is set
   * \param option_name name of the option
   * \param error_msg error message if option is missing
   *
   * CheckRequired is called from EvaluateOptions if a required options is set.
   * If not, the list of possible options is shown and an exception with
   * the error message given in error_msg is thrown
   */
  void CheckRequired(const std::string &option_name,
                     const std::string &error_msg = "");

  /**
   * \brief Main body of application
   *
   * Run is called after command line was parsed + evaluated. All
   * the work should be done in here.
   */
  virtual void Run() = 0;

  /**
   * \brief add option for command line
   * \param group group string
   * \return easy_init of boost, see documentation
   *
   * Adds an option to the available command line options. If no group is
   * specified, it is added to the standard group (Allowed Options). If group
   * is given, a sub group for this set of options will be created.
   */
  boost::program_options::options_description_easy_init AddProgramOptions(
      const std::string &group = "");

  /**
   * \brief get available program options & descriptions
   * \return variables_map (see boost documentation)
   */
  boost::program_options::variables_map &OptionsMap() { return _op_vm; }
  boost::program_options::options_description &OptionsDesc() {
    return _op_desc;
  }

  /**
   * \brief filters out the Hidden group from the options descriptions
   * \return Option descriptions without the "Hidden" group
   */
  boost::program_options::options_description &VisibleOptions() {
    return _visible_options;
  }

  /**
   * \brief call StopExecution after EvaluateOptions
   *
   * This is useful if the program executes an operation in EvaluateOptions
   * and then wants to stop execution successfully. Call StopExecution and
   * return true in EvaluateOptions.
   */
  void StopExecution() { _continue_execution = false; }

  /// length of the output help
  enum HelpType { HelpShort, HelpLong };

  /**
   * \brief Print long/short descriptions of calculators
   *
   * @param calculator_name name of a calculator
   * @param help_path path in VOTCASHARE were xml file with help is stored
   * @param helptype long or short (with options) help
   */
  void PrintDescription(std::ostream &out, const std::string &calculator_name,
                        const std::string help_path, HelpType helptype);

 protected:
  /// Variable map containing all program options
  boost::program_options::variables_map _op_vm;

  /// program options required by all applications
  boost::program_options::options_description _op_desc;

  std::map<std::string, boost::program_options::options_description> _op_groups;

  virtual void ShowHelpText(std::ostream &out);

  void ShowManPage(std::ostream &out);

  void ShowTEXPage(std::ostream &out);

  bool _continue_execution;

 private:
  /// get input parameters from file, location may be specified in command line
  void ParseCommandLine(int argc, char **argv);

  /// program options without the Hidden group
  boost::program_options::options_description _visible_options;
};

}  // namespace tools
}  // namespace votca

#endif /* __VOTCA_APPLICATION_H */
