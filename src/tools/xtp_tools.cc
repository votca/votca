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
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <votca/tools/property.h>
#include <votca/xtp/qmtool.h>
#include <votca/xtp/toolfactory.h>
#include <votca/xtp/version.h>
#include <votca/xtp/xtpapplication.h>

using namespace std;
using namespace votca;

class XtpTools : public xtp::XtpApplication {
 public:
  XtpTools() {}

  ~XtpTools() {}

  string ProgramName() { return "xtp_tools"; }

  void HelpText(ostream& out) {
    out << "Runs excitation/charge transport tools" << endl;
  }

  void SetTool(xtp::QMTool* tool) {
    _tool = std::unique_ptr<xtp::QMTool>(tool);
  }
  void Initialize();
  bool EvaluateOptions();
  void Run(void);

  void BeginEvaluate(int nThreads);
  bool Evaluate();

 private:
  tools::Property _options;
  std::unique_ptr<xtp::QMTool> _tool;
};

namespace propt = boost::program_options;

void XtpTools::Initialize() {

  xtp::QMToolFactory::RegisterAll();
  xtp::XtpApplication::Initialize();

  // Tools-related
  AddProgramOptions("Tools")("execute,e", propt::value<string>(),
                             "List of tools separated by ',' or ' '");
  AddProgramOptions("Tools")("list,l", "Lists all available tools");
  AddProgramOptions("Tools")("description,d", propt::value<string>(),
                             "Short description of a tool");
  // Options-related
  AddProgramOptions()("nthreads,t", propt::value<int>()->default_value(1),
                      "  number of threads to create");
}

bool XtpTools::EvaluateOptions() {

  string helpdir = "xtp/xml";

  if (OptionsMap().count("list")) {
    cout << "Available XTP tools: \n";

    for (const auto& tool : xtp::QMTools().getObjects()) {
      PrintDescription(std::cout, tool.first, helpdir, Application::HelpShort);
    }
    StopExecution();
    return true;
  }

  if (OptionsMap().count("description")) {
    CheckRequired("description", "no tool is given");
    tools::Tokenizer tok(OptionsMap()["description"].as<string>(), " ,\n\t");
    // loop over the names in the description string
    for (const std::string& n : tok) {
      // loop over tools
      bool printerror = true;
      for (const auto& tool : xtp::QMTools().getObjects()) {
        if (n.compare(tool.first) == 0) {
          PrintDescription(std::cout, tool.first, helpdir,
                           Application::HelpLong);
          printerror = false;
          break;
        }
      }
      if (printerror) cout << "Tool " << n << " does not exist\n";
    }
    StopExecution();
    return true;
  }

  CheckRequired("execute", "Nothing to do here: Abort.");
  CheckRequired("options", "Please provide an xml file with tool options");

  tools::Tokenizer xtools(OptionsMap()["execute"].as<string>(), " ,\n\t");
  std::vector<std::string> calc_string = xtools.ToVector();
  if (calc_string.size() != 1) {
    throw std::runtime_error(
        "You can only run one calculator at the same time.");
  }

  bool found_calc = false;
  for (const auto& tool : xtp::QMTools().getObjects()) {
    if (calc_string[0].compare(tool.first) == 0) {
      cout << " This is a XTP app" << endl;
      this->SetTool(xtp::QMTools().Create(calc_string[0]));
      found_calc = true;
      break;
    }
  }
  if (!found_calc) {
    cout << "Tool " << calc_string[0] << " does not exist\n";
    StopExecution();
  } else {
    cout << "Registered " << calc_string[0] << endl;
  }
  return 1;
}

void XtpTools::Run() {

  string optionsFile = _op_vm["options"].as<string>();
  tools::load_property_from_xml(_options, optionsFile);

  int nThreads = OptionsMap()["nthreads"].as<int>();
  std::string name = ProgramName();
  if (VersionString() != "") name = name + ", version " + VersionString();
  xtp::HelpTextHeader(name);
  cout << "Initializing tool " << endl;
  BeginEvaluate(nThreads);

  cout << "Evaluating tool " << endl;

  Evaluate();
}

void XtpTools::BeginEvaluate(int nThreads = 1) {
  cout << "... " << _tool->Identify() << " " << flush;
  _tool->setnThreads(nThreads);
  _tool->Initialize(_options);
  cout << endl;
}

bool XtpTools::Evaluate() {

  cout << "... " << _tool->Identify() << " " << flush;
  _tool->Evaluate();
  cout << endl;
  return true;
}

int main(int argc, char** argv) {

  XtpTools xtpapp;
  return xtpapp.Exec(argc, argv);
}
