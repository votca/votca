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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <votca/tools/property.h>
#include <votca/xtp/qmtool.h>
#include <votca/xtp/toolfactory.h>
#include <votca/xtp/version.h>
#include <votca/xtp/xtpapplication.h>

using namespace votca;

class XtpTools : public xtp::XtpApplication {
 public:
  XtpTools() = default;

  ~XtpTools() override = default;

  std::string ProgramName() override { return "xtp_tools"; }

  void HelpText(std::ostream& out) override {
    out << "Runs excitation/charge transport tools\n";
  }

  void SetTool(xtp::QMTool* tool) {
    _tool = std::unique_ptr<xtp::QMTool>(tool);
  }
  void Initialize() override;
  bool EvaluateOptions() override;
  void Run(void) override;

  void BeginEvaluate(Index nThreads);
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
  AddProgramOptions("Tools")("execute,e", propt::value<std::string>(),
                             "List of tools separated by ',' or ' '");
  AddProgramOptions("Tools")("list,l", "Lists all available tools");
  AddProgramOptions("Tools")("description,d", propt::value<std::string>(),
                             "Short description of a tool");
  AddProgramOptions("Tools")("name,n", propt::value<std::string>(),
                             "Name of the job to run");

  // Options-related
  AddProgramOptions()("nthreads,t", propt::value<Index>()->default_value(1),
                      "  number of threads to create");
}

bool XtpTools::EvaluateOptions() {

  std::string helpdir = "xtp/xml";

  if (OptionsMap().count("list")) {
    std::cout << "Available XTP tools: \n";

    for (const auto& tool : xtp::QMTools().getObjects()) {
      PrintDescription(std::cout, tool.first, helpdir, Application::HelpShort);
    }
    StopExecution();
    return true;
  }

  if (OptionsMap().count("description")) {
    CheckRequired("description", "no tool is given");
    tools::Tokenizer tok(OptionsMap()["description"].as<std::string>(),
                         " ,\n\t");
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
      if (printerror) {
        std::cout << "Tool " << n << " does not exist\n";
      }
    }
    StopExecution();
    return true;
  }

  CheckRequired("execute", "Please provide the name of the tool to execute");

  tools::Tokenizer xtools(OptionsMap()["execute"].as<std::string>(), " ,\n\t");
  std::vector<std::string> calc_string = xtools.ToVector();
  if (calc_string.size() != 1) {
    throw std::runtime_error(
        "You can only run one calculator at the same time.");
  }

  CheckRequired(
      "name", "Please provide the job name to run (same as the xyz file name)");

  bool found_calc = false;
  for (const auto& tool : xtp::QMTools().getObjects()) {
    if (calc_string[0].compare(tool.first) == 0) {
      this->SetTool(xtp::QMTools().Create(calc_string[0]));
      found_calc = true;
      break;
    }
  }
  if (!found_calc) {
    std::cout << "Tool " << calc_string[0] << " does not exist\n";
    StopExecution();
  } else {
    std::cout << "Registered " << calc_string[0];
  }
  return 1;
}

void XtpTools::Run() {

  auto it = _op_vm.find("options");
  if (it != _op_vm.cend()) {
    std::string optionsFile = _op_vm["options"].as<std::string>();
    _options.LoadFromXML(optionsFile);
  } else {
    // Empty user options
    tools::Property& opts = _options.add("options", "");
    opts.add(_tool->Identify(), "");
  }

  std::string job_name = _op_vm["name"].as<std::string>();
  tools::Property& opts = _options.get("options." + _tool->Identify());
  opts.add("job_name", job_name);

  Index nThreads = OptionsMap()["nthreads"].as<Index>();
  std::string name = ProgramName();
  if (VersionString() != "") {
    name = name + ", version " + VersionString();
  }
  xtp::HelpTextHeader(name);
  std::cout << "Initializing tool\n";
  BeginEvaluate(nThreads);

  std::cout << "Evaluating tool\n";

  Evaluate();
}

void XtpTools::BeginEvaluate(Index nThreads = 1) {
  std::cout << "... " << _tool->Identify() << " " << std::flush;
  _tool->setnThreads(nThreads);
  _tool->Initialize(_options);
}

bool XtpTools::Evaluate() {

  std::cout << "... " << _tool->Identify() << " " << std::flush;
  _tool->Evaluate();
  return true;
}

int main(int argc, char** argv) {

  XtpTools xtpapp;
  return xtpapp.Exec(argc, argv);
}
