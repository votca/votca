
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

// Local VOTCA includes
#include "votca/xtp/calculatorfactory.h"
#include "votca/xtp/stateapplication.h"
#include <memory>

using namespace votca;

class XtpRun final : public xtp::StateApplication {
 public:
  XtpRun() { xtp::Calculatorfactory::RegisterAll(); }

  ~XtpRun() = default;
  std::string ProgramName() final { return "xtp_run"; }

  void HelpText(std::ostream& out) final {
    out << "Runs excitation/charge transport calculators\n";
  }

 protected:
  void CreateCalculator(const std::string& name);
  void ConfigCalculator();
  bool savetoStateFile() const final { return _calc->WriteToStateFile(); }

  bool EvaluateFrame(votca::xtp::Topology& top) final;
  std::string CalculatorType() const { return "Calculator"; }
  void CheckOptions() final{};
  std::vector<std::string> CalculatorNames() const {
    return xtp::Calculators().getKeys();
  }

  void AddCommandLineOpt() final{};

 private:
  std::unique_ptr<xtp::QMCalculator> _calc = nullptr;
};

void XtpRun::CreateCalculator(const std::string& name) {
  _calc = xtp::Calculators().Create(name);
}

void XtpRun::ConfigCalculator() {
  std::cout << "... " << _calc->Identify() << std::endl;
  Index nThreads = OptionsMap()["nthreads"].as<Index>();
  _calc->setnThreads(nThreads);
  _calc->Initialize(_options);
}

bool XtpRun::EvaluateFrame(xtp::Topology& top) {
  std::cout << "... " << _calc->Identify() << std::endl;
  return _calc->EvaluateFrame(top);
}

int main(int argc, char** argv) {

  XtpRun xtprun;
  return xtprun.Exec(argc, argv);
}
