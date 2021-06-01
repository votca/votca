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

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/globals.h>
#include <votca/tools/propertyiomanipulator.h>

// Local VOTCA includes
#include "votca/xtp/version.h"
#include "votca/xtp/xtpapplication.h"
#include "votca/xtp/openmp_cuda.h"

namespace votca {
namespace xtp {
XtpApplication::XtpApplication() = default;
XtpApplication::~XtpApplication() = default;

/**
 * \brief Adds program options to the executable
 *
 * Every executable requires option file for calculators it is running
 * It is thus a part of the base XtpApplication class
 *
 */
void XtpApplication::Initialize(void) {

  namespace propt = boost::program_options;
  AddProgramOptions()(
      "options,o", propt::value<std::string>(),
      std::string("  " + CalculatorType() + " options").c_str());

  AddProgramOptions()("nthreads,t", propt::value<Index>()->default_value(1),
                      "  number of threads to create");

  AddProgramOptions()(
      "execute,e", propt::value<std::string>(),
      std::string("Name of " + CalculatorType() + " to run").c_str());
  AddProgramOptions()(
      "list,l",
      std::string("Lists all available " + CalculatorType() + "s").c_str());
  AddProgramOptions()(
      "description,d", propt::value<std::string>(),
      std::string("Short description of a " + CalculatorType() + "s").c_str());

#ifdef USE_CUDA
  AddProgramOptions()("gpus,g", propt::value<Index>()->default_value(-1),
                      "  Number of gpus to use");
#endif

  AddCommandLineOptions();
}

bool XtpApplication::EvaluateOptions() {

  if (OptionsMap().count("list")) {
    std::cout << "Available XTP" + CalculatorType() + "s: \n";
    for (const auto& name : CalculatorNames()) {
      PrintDescription(std::cout, name, "xtp/xml", Application::HelpShort);
    }
    StopExecution();
    return true;
  }
#ifdef USE_CUDA
  OpenMP_CUDA::SetNoGPUs(OptionsMap()["gpus"].as<Index>());
#endif

  if (OptionsMap().count("description")) {
    CheckRequired("description", "no " + CalculatorType() + " is given");
    tools::Tokenizer tok(OptionsMap()["description"].as<std::string>(),
                         " ,\n\t");
    for (const std::string& n : tok) {
      if (CalcExists(n)) {
        PrintDescription(std::cout, n, "xtp/xml", Application::HelpLong);
      } else {
        std::cout << CalculatorType() << " " << n << " does not exist\n";
      }
    }
    StopExecution();
    return true;
  }
  CheckRequired("execute", "Nothing to do here: Abort.");

  tools::Tokenizer calcs(OptionsMap()["execute"].as<std::string>(), " ,\n\t");
  std::vector<std::string> calc_string = calcs.ToVector();
  if (calc_string.size() != 1) {
    throw std::runtime_error("You can only run one " + CalculatorType() +
                             " at the same time.");
  }

  if (CalcExists(calc_string[0])) {
    CreateCalculator(calc_string[0]);
  } else {
    std::cout << CalculatorType() << " " << calc_string[0]
              << " does not exist\n";
    StopExecution();
  }

  EvaluateSpecificOptions();

  return true;
}

void XtpApplication::Run() {

  std::string name = ProgramName();
  if (VersionString() != "") {
    name = name + ", version " + VersionString();
  }
  xtp::HelpTextHeader(name);

  execute();
}

void XtpApplication::ShowHelpText(std::ostream& out) {
  std::string name = ProgramName();
  if (VersionString() != "") {
    name = name + ", version " + VersionString();
  }
  xtp::HelpTextHeader(name);
  HelpText(out);
  out << "\n\n" << VisibleOptions() << std::endl;
}

}  // namespace xtp
}  // namespace votca
