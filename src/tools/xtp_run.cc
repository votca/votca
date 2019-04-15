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
#include <votca/ctp/calculatorfactory.h>
#include <votca/xtp/calculatorfactory.h>
#include <votca/xtp/sqlapplication.h>

using namespace std;
using namespace votca;

class XtpRun : public xtp::SqlApplication {
 public:
  string ProgramName() { return "xtp_run"; }

  void HelpText(ostream& out) {
    out << "Runs excitation/charge transport calculators" << endl;
  }

  void HelpText(){};

  void Initialize();
  bool EvaluateOptions();

 private:
  // void    PrintDescription(string name, HelpOutputType _help_output_type);
};

namespace propt = boost::program_options;

void XtpRun::Initialize() {
  xtp::Calculatorfactory::RegisterAll();
  ctp::Calculatorfactory::RegisterAll();
  xtp::SqlApplication::Initialize();

  AddProgramOptions("Calculators")(
      "execute,e", propt::value<string>(),
      "List of calculators separated by ',' or ' '");
  AddProgramOptions("Calculators")("list,l", "Lists all available calculators");
  AddProgramOptions("Calculators")("description,d", propt::value<string>(),
                                   "Short description of a calculator");
}

bool XtpRun::EvaluateOptions() {

  string helpdir = "xtp/xml";
  string ctphelpdir = "ctp/xml";
  if (OptionsMap().count("list")) {
    cout << "Available XTP calculators: \n";
    for (const auto& calc : xtp::Calculators().getObjects()) {
      PrintDescription(std::cout, calc.first, helpdir, Application::HelpShort);
    }
    cout << "Available (wrapped) CTP calculators: \n";
    for (const auto& calc : ctp::Calculators().getObjects()) {
      bool printctp = true;
      std::string ctpcalc = calc.first.c_str();
      for (const auto& xtpcalc : xtp::Calculators().getObjects()) {
        if (ctpcalc.compare(xtpcalc.first.c_str()) == 0) {
          printctp = false;
          break;
        }
      }
      if (printctp) {
        PrintDescription(std::cout, calc.first, ctphelpdir,
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
      for (const auto& calc : xtp::Calculators().getObjects()) {

        if (n.compare(calc.first.c_str()) == 0) {
          PrintDescription(std::cout, calc.first, helpdir,
                           Application::HelpLong);
          printerror = false;
          break;
        }
      }
      for (const auto& calc : ctp::Calculators().getObjects()) {

        if (n.compare(calc.first.c_str()) == 0) {
          bool printctp = true;
          std::string ctpcalc = calc.first.c_str();
          for (const auto& xtpcalc : xtp::Calculators().getObjects()) {
            if (ctpcalc.compare(xtpcalc.first.c_str()) == 0) {
              printctp = false;
              break;
            }
          }
          if (printctp) {
            PrintDescription(std::cout, calc.first, "ctp/xml",
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

  xtp::SqlApplication::EvaluateOptions();
  CheckRequired("options",
                "Please provide an xml file with calculator options");
  CheckRequired("execute", "Nothing to do here: Abort.");

  Tokenizer calcs(OptionsMap()["execute"].as<string>(), " ,\n\t");
  for (const std::string& n : calcs) {
    bool _found_calc = false;
    for (const auto& calc : xtp::Calculators().getObjects()) {

      if (n.compare(calc.first.c_str()) == 0) {
        cout << " This is a XTP app" << endl;
        xtp::SqlApplication::AddCalculator(
            xtp::Calculators().Create(n.c_str()));
        _found_calc = true;
      }
    }

    if (!_found_calc) {
      for (const auto& calc : ctp::Calculators().getObjects()) {

        if (n.compare(calc.first.c_str()) == 0) {
          _found_calc = true;
          cout << " This is a CTP app" << endl;
          xtp::SqlApplication::AddCalculator(
              ctp::Calculators().Create(n.c_str()));
        }
      }
    }
    if (!_found_calc) {
      cout << "Calculator " << n << " does not exist\n";
      StopExecution();
    } else {
      load_property_from_xml(_options, _op_vm["options"].as<string>());
    }
  }
  return true;
}

int main(int argc, char** argv) {

  XtpRun xtprun;
  return xtprun.Exec(argc, argv);
}
