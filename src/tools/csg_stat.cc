/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#include "csg_stat_imc.h"
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <votca/csg/csgapplication.h>
#include <votca/csg/version.h>

// using namespace votca::tools;
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

  CsgApplication::Worker *ForkWorker() override { return _imc.ForkWorker(); }

  void MergeWorker(CsgApplication::Worker *worker) override {
    _imc.MergeWorker(worker);
  }

 public:
  Imc _imc;
  int _block_length;
  string _extension;
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
      "block-length", boost::program_options::value<int>(),
      "  write blocks of this length, the averages are cleared after every "
      "write")("ext",
               boost::program_options::value<string>(&_extension)
                   ->default_value("dist.new"),
               "Extension of the output");
}

bool CsgStatApp::EvaluateOptions() {
  CsgApplication::EvaluateOptions();
  CheckRequired("options");
  CheckRequired("trj", "no trajectory file specified");

  _imc.LoadOptions(OptionsMap()["options"].as<string>());

  if (OptionsMap().count("block-length")) {
    _imc.BlockLength(OptionsMap()["block-length"].as<int>());
  } else {
    _imc.BlockLength(0);
  }

  if (OptionsMap().count("do-imc")) {
    _imc.DoImc(true);
  }

  _imc.Extension(_extension);

  _imc.Initialize();
  return true;
}

void CsgStatApp::BeginEvaluate(Topology *top, Topology *top_ref) {
  _imc.BeginEvaluate(top, top_ref);
}

void CsgStatApp::EndEvaluate() { _imc.EndEvaluate(); }

int main(int argc, char **argv) {
  CsgStatApp app;
  return app.Exec(argc, argv);
}
