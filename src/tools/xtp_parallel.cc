/*
 *            Copyright 2009-2017 The VOTCA Development Team
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

#include <votca/ctp/jobcalculatorfactory.h>

using namespace std;
using namespace votca;

class XtpParallel : public xtp::JobApplication {
 public:
  string ProgramName() { return "xtp_parallel"; }

  void HelpText(ostream& out) {
    out << "Runs job-based heavy-duty calculators" << endl;
  }
  void HelpText(){};

  void Initialize();
  bool EvaluateOptions();

 private:
  // void    PrintDescription(string name, HelpOutputType _help_output_type);
};

namespace propt = boost::program_options;

void XtpParallel::Initialize() {
  xtp::JobCalculatorfactory::RegisterAll();
  ctp::JobCalculatorfactory::RegisterAll();
  xtp::JobApplication::Initialize();

  AddProgramOptions("Calculators")(
      "execute,e", propt::value<string>(),
      "List of calculators separated by ',' or ' '");
  AddProgramOptions("Calculators")("list,l", "Lists all available calculators");
  AddProgramOptions("Calculators")("description,d", propt::value<string>(),
                                   "Short description of a calculator");
}

bool XtpParallel::EvaluateOptions() {

  if (OptionsMap().count("list")) {
    cout << "Available XTP calculators: \n";
    for (const auto& jobcalc : xtp::JobCalculators().getObjects()) {
      PrintDescription(std::cout, jobcalc.first, "xtp/xml",
                       Application::HelpShort);
    }

    cout << "Available (wrapped) CTP calculators: \n";
    for (const auto& jobcalc : ctp::JobCalculators().getObjects()) {
      bool printctp = true;
      std::string ctpcalc = jobcalc.first.c_str();
      for (const auto& xtpjobcalc : xtp::JobCalculators().getObjects()) {
        if (ctpcalc.compare(xtpjobcalc.first.c_str()) == 0) {
          printctp = false;
          break;
        }
      }

      if (printctp) {
        PrintDescription(std::cout, jobcalc.first, "ctp/xml",
                         Application::HelpShort);
      }
    }
    StopExecution();
    return true;
  }

  if (OptionsMap().count("description")) {
    CheckRequired("description", "no calculator is given");
    Tokenizer tok(OptionsMap()["description"].as<string>(), " ,\n\t");
    // loop over the names in the description string
    for (const std::string& n : tok) {
      // loop over calculators
      bool printerror = true;
      for (const auto& jobcalc : xtp::JobCalculators().getObjects()) {

        if (n.compare(jobcalc.first.c_str()) == 0) {
          PrintDescription(std::cout, jobcalc.first, "xtp/xml",
                           Application::HelpLong);
          printerror = false;
          break;
        }
      }
      for (const auto& jobcalc : ctp::JobCalculators().getObjects()) {

        if (n.compare(jobcalc.first.c_str()) == 0) {
          bool printctp = true;
          std::string ctpcalc = jobcalc.first.c_str();
          for (const auto& xtpjobcalc : xtp::JobCalculators().getObjects()) {
            if (ctpcalc.compare(xtpjobcalc.first.c_str()) == 0) {
              printctp = false;
              break;
            }
          }
          if (printctp) {
            PrintDescription(std::cout, jobcalc.first, "ctp/xml",
                             Application::HelpLong);
            printerror = false;
            break;
          }
        }
      }
      if (printerror) cout << "Calculator " << n << " does not exist\n";
    }
    StopExecution();
    return true;
  }

  xtp::JobApplication::EvaluateOptions();
  CheckRequired("execute", "Nothing to do here: Abort.");

  Tokenizer calcs(OptionsMap()["execute"].as<string>(), " ,\n\t");
  for (const std::string& n : calcs) {
    bool _found_calc = false;
    for (const auto& jobcalc : xtp::JobCalculators().getObjects()) {

      if (n.compare(jobcalc.first.c_str()) == 0) {
        cout << " This is a XTP app" << endl;
        xtp::JobApplication::AddCalculator(
            xtp::JobCalculators().Create(n.c_str()));
        _found_calc = true;
      }
    }
    if (!_found_calc) {
      for (const auto& jobcalc : ctp::JobCalculators().getObjects()) {

        if (n.compare(jobcalc.first.c_str()) == 0) {
          cout << " This is a CTP app" << endl;
          xtp::JobApplication::AddCalculator(
              ctp::JobCalculators().Create(n.c_str()));
          _found_calc = true;
        }
      }
    }
    if (!_found_calc) {
      cout << "Jobcalculator " << n << " does not exist\n";
      StopExecution();
    }
  }

  return true;
}

int main(int argc, char** argv) {

  XtpParallel xtprun;
  return xtprun.Exec(argc, argv);
}
