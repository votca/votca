/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Standard includes
#include <cstdlib>
#include <fstream>
#include <memory>

// Third party includes
#include <boost/program_options.hpp>

// Local VOTCA includes
#include "votca/csg/csgapplication.h"
#include "votca/csg/version.h"

// Local private VOTCA includes
#include "csg_stat_imc.h"

using namespace std;
using namespace votca::csg;

class CsgStatApp : public CsgApplication {
 public:
  string ProgramName() override { return "csg_stat"; }
  void HelpText(ostream &out) override;

  bool DoTrajectory() override { return true; }
  bool DoMapping() override { return true; }
  bool DoMappingDefault(void) override { return false; }
  bool DoThreaded() override { return true; }
  bool SynchronizeThreads() override { return true; }
  void Initialize() override;
  bool EvaluateOptions() override;

  void BeginEvaluate(Topology *top, Topology *top_ref) override;
  void EndEvaluate() override;

  std::unique_ptr<CsgApplication::Worker> ForkWorker() override {
    return imc_.ForkWorker();
  }

  void MergeWorker(CsgApplication::Worker *worker) override {
    imc_.MergeWorker(worker);
  }

 public:
  Imc imc_;
  votca::Index block_length_;
  string extension_;
};

void CsgStatApp::HelpText(ostream &out) {
  out << "Calculate all distributions (bonded and non-bonded) specified in "
         "options file.\n"
         "Optionally calculates update Eigen::Matrix3d for invere Monte Carlo. "
         "This "
         "program\n"
         "is called inside the inverse scripts. Unlike csg_boltzmann, big "
         "systems\n"
         "can be treated as well as non-bonded interactions can be evaluated.";
}

void CsgStatApp::Initialize() {
  CsgApplication::Initialize();
  AddProgramOptions("Specific options")("options",
                                        boost::program_options::value<string>(),
                                        "  options file for coarse graining")(
      "do-imc", "  write out additional Inverse Monte Carlo data")(
      "only-intra-nb", "  only do intramolecular non-bonded interactions")(
      "block-length", boost::program_options::value<votca::Index>(),
      "  write blocks of this length, the averages are cleared after every "
      "write")("ext",
               boost::program_options::value<string>(&extension_)
                   ->default_value("dist.new"),
               "Extension of the output");
}

bool CsgStatApp::EvaluateOptions() {
  CsgApplication::EvaluateOptions();
  CheckRequired("options");
  CheckRequired("trj", "no trajectory file specified");

  imc_.LoadOptions(OptionsMap()["options"].as<string>());

  if (OptionsMap().count("block-length")) {
    imc_.BlockLength(OptionsMap()["block-length"].as<votca::Index>());
  } else {
    imc_.BlockLength(0);
  }

  if (OptionsMap().count("do-imc")) {
    imc_.DoImc(true);
  }

  if (OptionsMap().count("only-intra-nb")) {
    imc_.OnlyIntraNB(true);
  }

  imc_.Extension(extension_);

  imc_.Initialize();
  return true;
}

void CsgStatApp::BeginEvaluate(Topology *top, Topology *top_ref) {
  imc_.BeginEvaluate(top, top_ref);
}

void CsgStatApp::EndEvaluate() { imc_.EndEvaluate(); }

int main(int argc, char **argv) {
  CsgStatApp app;
  return app.Exec(argc, argv);
}
