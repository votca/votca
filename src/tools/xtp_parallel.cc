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
  bool savetoStateFile() const final { return import_; }

  bool EvaluateFrame(votca::xtp::Topology& top) final;
  std::string CalculatorType() const { return "Calculator"; }
  void CheckOptions() final;
  std::vector<std::string> CalculatorNames() const {
    return xtp::JobCalculators().getKeys();
  }

  void AddCommandLineOpt() final;

 private:
  bool generate_input_ = false;
  bool run_ = false;
  bool import_ = false;
  xtp::ProgObserver<std::vector<xtp::Job>> progObs_;
  std::unique_ptr<xtp::JobCalculator> calc_ = nullptr;
};

void XtpParallel::AddCommandLineOpt() {
  namespace propt = boost::program_options;
  AddProgramOptions()("ompthreads,x", propt::value<Index>()->default_value(1),
                      "  number of openmp threads to create in each thread");
  AddProgramOptions()("restart,r",
                      propt::value<std::string>()->default_value(""),
                      "  restart pattern: 'host(pc1:234) stat(FAILED)'");
  AddProgramOptions()("cache,q", propt::value<Index>()->default_value(8),
                      "  assigns jobs in blocks of this size");
  AddProgramOptions()("jobs,j",
                      propt::value<std::string>()->default_value("run"),
                      "  task(s) to perform: write, run, read");
  AddProgramOptions()("maxjobs,m", propt::value<Index>()->default_value(-1),
                      "  maximum number of jobs to process (-1 = inf)");
}

void XtpParallel::CheckOptions() {
  std::string jobstr = OptionsMap()["jobs"].as<std::string>();
  generate_input_ = (jobstr == "write");
  run_ = (jobstr == "run");
  import_ = (jobstr == "read");
}

void XtpParallel::CreateCalculator(const std::string& name) {
  calc_ = xtp::JobCalculators().Create(name);
}

void XtpParallel::ConfigCalculator() {
  std::cout << "... " << calc_->Identify() << std::endl;

  progObs_.InitCmdLineOpts(OptionsMap());
  calc_->setnThreads(OptionsMap()["nthreads"].as<Index>());
  calc_->setOpenMPThreads(OptionsMap()["ompthreads"].as<Index>());
  calc_->setProgObserver(&progObs_);
  calc_->Initialize(options_);
  std::cout << std::endl;
}

bool XtpParallel::EvaluateFrame(xtp::Topology& top) {
  std::cout << "... " << calc_->Identify() << " " << std::flush;
  if (generate_input_) {
    calc_->WriteJobFile(top);
  } else if (run_) {
    calc_->EvaluateFrame(top);
  } else if (import_) {
    calc_->ReadJobFile(top);
  }
  std::cout << std::endl;
  return true;
}

int main(int argc, char** argv) {

  XtpParallel xtprun;
  return xtprun.Exec(argc, argv);
}
