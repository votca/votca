/*
 *            Copyright 2009-2021 The VOTCA Development Team
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

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/qmtool.h"
#include "votca/xtp/toolfactory.h"
#include "votca/xtp/version.h"
#include "votca/xtp/xtpapplication.h"

using namespace votca;

class XtpTools final : public xtp::XtpApplication {
 public:
  XtpTools(){xtp::QMToolFactory::RegisterAll(); }

  ~XtpTools() = default;

  std::string ProgramName() final { return "xtp_tools"; }

  void HelpText(std::ostream& out) final {
    out << "Runs excitation/charge transport tools\n";
  }

protected:

  void CreateCalculator(const std::string& name) final;

  void execute() final;
  std::string CalculatorType() const final{return "Tool";}
  void EvaluateSpecificOptions() final;
  std::vector<std::string> CalculatorNames() const final{return xtp::QMTools().getKeys();}

  void AddCommandLineOptions() final;

 private:
  std::unique_ptr<xtp::QMTool> _tool;
};

void XtpTools::CreateCalculator(const std::string& name){
_tool=xtp::QMTools().Create(name);
}
void XtpTools::AddCommandLineOptions() {
  namespace propt = boost::program_options;
  AddProgramOptions()("name,n", propt::value<std::string>(),
                      "Name of the job to run");
}

void XtpTools::EvaluateSpecificOptions() {
  CheckRequired(
      "name", "Please provide the job name to run (same as the xyz file name)");
}

void XtpTools::execute() {

  if (_op_vm.find("options") != _op_vm.cend()) {
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
  
  std::cout << "Initializing tool\n";
    std::cout << "... " << _tool->Identify() << " " << std::flush;
  _tool->setnThreads(nThreads);
  _tool->Initialize(_options);

  std::cout << "Evaluating tool\n";
  std::cout << "... " << _tool->Identify() << " " << std::flush;
  _tool->Evaluate();
}


int main(int argc, char** argv) {

  XtpTools xtpapp;
  return xtpapp.Exec(argc, argv);
}
