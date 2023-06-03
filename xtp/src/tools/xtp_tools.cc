/*
 *            Copyright 2009-2023 The VOTCA Development Team
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
  XtpTools() { xtp::QMToolFactory(); }

  ~XtpTools() = default;

  std::string ProgramName() final { return "xtp_tools"; }

  void HelpText(std::ostream& out) final {
    out << "Runs excitation/charge transport tools\n";
  }

 protected:
  void CreateCalculator(const std::string& name) final;

  void execute() final;
  std::string CalculatorType() const final { return "Tool"; }
  void EvaluateSpecificOptions() final;
  std::vector<std::string> CalculatorNames() const final {
    return xtp::QMToolFactory().getKeys();
  }

  void AddCommandLineOptions() final;

 private:
  std::unique_ptr<xtp::QMTool> tool_;
};

void XtpTools::CreateCalculator(const std::string& name) {
  tool_ = xtp::QMToolFactory().Create(name);
}
void XtpTools::AddCommandLineOptions() {}

void XtpTools::EvaluateSpecificOptions() {}

void XtpTools::execute() {

  Index nThreads = OptionsMap()["nthreads"].as<Index>();

  std::cout << "Initializing tool\n";
  std::cout << "... " << tool_->Identify() << " " << std::flush;
  tool_->setnThreads(nThreads);
  tool_->Initialize(options_);

  std::cout << "Evaluating tool\n";
  std::cout << "... " << tool_->Identify() << " " << std::flush;
  tool_->Evaluate();
}

int main(int argc, char** argv) {

  XtpTools xtpapp;
  return xtpapp.Exec(argc, argv);
}
