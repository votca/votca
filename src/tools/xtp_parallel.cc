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

#include <iostream>
#include <stdlib.h>
#include <string>
#include <votca/xtp/jobapplication.h>
#include <votca/xtp/jobcalculatorfactory.h>

using namespace std;
using namespace votca;

class XtpParallel : public xtp::JobApplication {
 public:
  string ProgramName() override { return "xtp_parallel"; }

  void HelpText(ostream& out) override {
    out << "Runs job-based heavy-duty calculators" << endl;
  }
  void HelpText(){};

  void Initialize() override;
  bool EvaluateOptions() override;

 private:
  // void    PrintDescription(string name, HelpOutputType _help_output_type);
};

namespace propt = boost::program_options;

void XtpParallel::Initialize() {
  xtp::JobCalculatorfactory::RegisterAll();
  xtp::JobApplication::Initialize();

  AddProgramOptions("Calculator")("execute,e", propt::value<string>(),
                                  "Name of calculator to run");
  AddProgramOptions("Calculator")("list,l", "Lists all available calculators");
  AddProgramOptions("Calculator")("description,d", propt::value<string>(),
                                  "Short description of a calculator");
}

bool XtpParallel::EvaluateOptions() {

  if (OptionsMap().count("list")) {
    cout << "Available XTP calculators: \n";
    for (const auto& jobcalc : xtp::JobCalculators().getObjects()) {
      PrintDescription(std::cout, jobcalc.first, "xtp/xml",
                       Application::HelpShort);
    }
    StopExecution();
    return true;
  }

  if (OptionsMap().count("description")) {
    CheckRequired("description", "no calculator is given");
    tools::Tokenizer tok(OptionsMap()["description"].as<string>(), " ,\n\t");
    // loop over the names in the description string
    for (const std::string& n : tok) {
      // loop over calculators
      bool printerror = true;
      for (const auto& jobcalc : xtp::JobCalculators().getObjects()) {
        if (n.compare(jobcalc.first) == 0) {
          PrintDescription(std::cout, jobcalc.first, "xtp/xml",
                           Application::HelpLong);
          printerror = false;
          break;
        }
      }
      if (printerror) cout << "Calculator " << n << " does not exist\n";
    }
    StopExecution();
    return true;
  }

  xtp::JobApplication::EvaluateOptions();
  CheckRequired("execute", "Nothing to do here: Abort.");

  tools::Tokenizer calcs(OptionsMap()["execute"].as<string>(), " ,\n\t");
  std::vector<std::string> calc_string = calcs.ToVector();
  if (calc_string.size() != 1) {
    throw std::runtime_error(
        "You can only run one calculator at the same time.");
  }

  bool found_calc = false;
  for (const auto& jobcalc : xtp::JobCalculators().getObjects()) {
    if (calc_string[0].compare(jobcalc.first) == 0) {
      cout << " This is a XTP app" << endl;
      xtp::JobApplication::SetCalculator(
          xtp::JobCalculators().Create(calc_string[0]));
      found_calc = true;
      break;
    }
  }

  if (!found_calc) {
    cout << "Jobcalculator " << calc_string[0] << " does not exist\n";
    StopExecution();
  }

  return true;
}

int main(int argc, char** argv) {

  XtpParallel xtprun;
  return xtprun.Exec(argc, argv);
}
