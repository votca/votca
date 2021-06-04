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

// Local VOTCA includes
#include "votca/xtp/jobcalculator.h"
#include "votca/xtp/jobcalculatorfactory.h"
#include "votca/xtp/progressobserver.h"
#include "votca/xtp/stateapplication.h"

using namespace votca;

class XtpParallel final : public xtp::StateApplication {
 public:
  XtpParallel() { xtp::JobCalculatorfactory::RegisterAll(); }

  ~XtpParallel() = default;
  std::string ProgramName() final { return "xtp_parallel"; }

  void HelpText(std::ostream& out) final {
    out << "Runs job-based heavy-duty calculators\n";
  }

 protected:
  void CreateCalculator(const std::string& name);
  void ConfigCalculator();
  bool savetoStateFile() const final { return _import; }

  bool EvaluateFrame(votca::xtp::Topology& top) final;
  std::string CalculatorType() const { return "Calculator"; }
  void CheckOptions() final;
  std::vector<std::string> CalculatorNames() const {
    return xtp::JobCalculators().getKeys();
  }

  void AddCommandLineOpt() final;

 private:
  bool _generate_input = false;
  bool _run = false;
  bool _import = false;
  xtp::ProgObserver<std::vector<xtp::Job>> _progObs;
  std::unique_ptr<xtp::JobCalculator> _calc = nullptr;
};

void XtpParallel::AddCommandLineOpt() {
  namespace propt = boost::program_options;
  AddProgramOptions()("ompthreads,p", propt::value<Index>()->default_value(1),
                      "  number of openmp threads to create in each thread");
  AddProgramOptions()("restart,r",
                      propt::value<std::string>()->default_value(""),
                      "  restart pattern: 'host(pc1:234) stat(FAILED)'");
  AddProgramOptions()("cache,c", propt::value<Index>()->default_value(8),
                      "  assigns jobs in blocks of this size");
  AddProgramOptions()("jobs,j",
                      propt::value<std::string>()->default_value("run"),
                      "  task(s) to perform: write, run, read");
  AddProgramOptions()("maxjobs,m", propt::value<Index>()->default_value(-1),
                      "  maximum number of jobs to process (-1 = inf)");
}

void XtpParallel::CheckOptions() {
  std::string jobstr = OptionsMap()["jobs"].as<std::string>();
  _generate_input = (jobstr == "write");
  _run = (jobstr == "run");
  _import = (jobstr == "read");
}

void XtpParallel::CreateCalculator(const std::string& name) {
  _calc = xtp::JobCalculators().Create(name);
}

void XtpParallel::ConfigCalculator() {
  std::cout << "... " << _calc->Identify() << std::endl;

  _progObs.InitCmdLineOpts(OptionsMap());
  _calc->setnThreads(OptionsMap()["nthreads"].as<Index>());
  _calc->setOpenMPThreads(OptionsMap()["ompthreads"].as<Index>());
  _calc->setProgObserver(&_progObs);
  _calc->Initialize(_options);
  std::cout << std::endl;
}

bool XtpParallel::EvaluateFrame(xtp::Topology& top) {
  std::cout << "... " << _calc->Identify() << " " << std::flush;
  if (_generate_input) {
    _calc->WriteJobFile(top);
  } else if (_run) {
    _calc->EvaluateFrame(top);
  } else if (_import) {
    _calc->ReadJobFile(top);
  }
  std::cout << std::endl;
  return true;
}

int main(int argc, char** argv) {

  XtpParallel xtprun;
  return xtprun.Exec(argc, argv);
}
