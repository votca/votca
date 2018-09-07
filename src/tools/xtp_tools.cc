/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/tools/property.h>
#include <votca/ctp/toolfactory.h>
#include <votca/xtp/toolfactory.h>
#include <votca/xtp/xtpapplication.h>
#include <stdio.h>
#include <votca/xtp/version.h>



using namespace std;
using namespace votca;

class XtpTools : public xtp::XtpApplication {
public:

  XtpTools() {
  }
  
  ~XtpTools(){
    for (auto* tool:_tools){
      delete tool;
    }
  }

  string ProgramName() {
    return "xtp_tools";
  }

  void HelpText(ostream &out) {
    out << "Runs excitation/charge transport tools" << endl;
  }

  void AddTool(votca::ctp::QMTool *tool) {
    _tools.push_back(tool);
  }
  void Initialize();
  bool EvaluateOptions();
  void Run(void);

  void BeginEvaluate(int nThreads);
  bool Evaluate();
  void EndEvaluate();

private:

  tools::Property _options;
  list< ctp::QMTool* > _tools;

};

namespace propt = boost::program_options;

void XtpTools::Initialize() {

  xtp::QMToolFactory::RegisterAll();
  ctp::QMToolFactory::RegisterAll();
  xtp::XtpApplication::Initialize();

  // Tools-related
  AddProgramOptions("Tools") ("execute,e", propt::value<string>(),
          "List of tools separated by ',' or ' '");
  AddProgramOptions("Tools") ("list,l",
          "Lists all available tools");
  AddProgramOptions("Tools") ("description,d", propt::value<string>(),
          "Short description of a tool");
  // Options-related
  AddProgramOptions() ("nthreads,t", propt::value<int>()->default_value(1),
          "  number of threads to create");

}

bool XtpTools::EvaluateOptions() {

  string helpdir = "xtp/xml";
  string ctphelpdir = "ctp/xml";

  if (OptionsMap().count("list")) {
    cout << "Available XTP tools: \n";
    for (const auto& tool:xtp::QMTools().getObjects()) {
      PrintDescription(std::cout, tool.first, helpdir, Application::HelpShort);
    }
    cout << "Available CTP tools: \n";
    for (const auto& tool:ctp::QMTools().getObjects()) {
      bool printctp = true;
      std::string ctpcalc = (tool.first).c_str();
      for (const auto& xtptool:xtp::QMTools().getObjects()) {
        if (ctpcalc.compare((xtptool.first).c_str()) == 0) {
          printctp = false;
          break;
        }
      }
      if (printctp) {
        PrintDescription(std::cout, tool.first, ctphelpdir, Application::HelpShort);
      }
    }
    StopExecution();
    return true;
  }


  if (OptionsMap().count("description")) {
    CheckRequired("description", "no tool is given");
    tools::Tokenizer tok(OptionsMap()["description"].as<string>(), " ,\n\t");
    // loop over the names in the description string
    for (const std::string& n :tok) {
      // loop over tools
      bool printerror = true;
      for ( const auto& tool:xtp::QMTools().getObjects()) {
        if (n.compare(tool.first.c_str()) == 0) {
          PrintDescription(std::cout, tool.first, helpdir, Application::HelpLong);
          printerror = false;
          break;
        }
      }
      for (const auto& tool:ctp::QMTools().getObjects()) {
        if (n.compare(tool.first.c_str()) == 0) {
          bool printctp = true;
          std::string ctpcalc = tool.first.c_str();
          for (const auto& xtptool:xtp::QMTools().getObjects()) {
            if (ctpcalc.compare(xtptool.first.c_str()) == 0) {
              printctp = false;
              break;
            }
          }
          if (printctp) {
            PrintDescription(std::cout, tool.first, "ctp/xml", Application::HelpLong);
            printerror = false;
            break;
          }
        }
      }
      if (printerror) cout << "Tool " << n << " does not exist\n";
    }
    StopExecution();
    return true;
  }
  
  Application::EvaluateOptions();
  CheckRequired("execute", "Nothing to do here: Abort.");
  CheckRequired("options", "Please provide an xml file with tool options");

  tools::Tokenizer xtools(OptionsMap()["execute"].as<string>(), " ,\n\t");
  for (const std::string& n :xtools) {
    bool _found_calc = false;
    for (const auto& tool:xtp::QMTools().getObjects()) {
      if (n.compare(tool.first.c_str()) == 0) {
        cout << " This is a XTP app" << endl;
        this->AddTool(xtp::QMTools().Create(n.c_str()));
        _found_calc = true;
      }
    }
    if (!_found_calc) {
      for (const auto& tool:ctp::QMTools().getObjects()) {

        if (n.compare(tool.first.c_str()) == 0) {
          cout << " This is a CTP app" << endl;
           this->AddTool(ctp::QMTools().Create(n.c_str()));
            _found_calc = true;
        }
      }
    }
    if (!_found_calc) {
      cout << "Tool " << n << " does not exist\n";
      StopExecution();
    }else{
      cout << "Registered " << n << endl;
    }
  }
  return 1;
}

void XtpTools::Run() {

  string optionsFile = _op_vm["options"].as<string>();
  tools::load_property_from_xml(_options, optionsFile);

  int nThreads = OptionsMap()["nthreads"].as<int>();
  std::string name = ProgramName();
  if (VersionString() != "") name = name + ", version " + VersionString();
  votca::xtp::HelpTextHeader(name);
  cout << "Initializing tools " << endl;
  BeginEvaluate(nThreads);

  cout << "Evaluating tools " << endl;
  
  Evaluate();
  EndEvaluate();
}

void XtpTools::BeginEvaluate(int nThreads = 1) {
  for (ctp::QMTool* tool: _tools) {
    cout << "... " << tool->Identify() << " " << flush;
    tool->setnThreads(nThreads);
    tool->Initialize(&_options);
    cout << endl;
  }
}

bool XtpTools::Evaluate() {
for (ctp::QMTool* tool: _tools) {
    cout << "... " << tool->Identify() << " " << flush;
    tool->Evaluate();
    cout << endl;
  }

  return true;
}

void XtpTools::EndEvaluate() {
for (ctp::QMTool* tool: _tools) {
    tool->EndEvaluate();
  }
}

int main(int argc, char** argv) {

  XtpTools xtpapp;
  return xtpapp.Exec(argc, argv);

}
