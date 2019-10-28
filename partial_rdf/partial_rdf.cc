/*
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#include "rdf_calculator.h"
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <votca/csg/cgengine.h>
#include <votca/csg/version.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

#include <stdlib.h>
#include <votca/csg/csgapplication.h>

// using namespace votca::tools;
using namespace std;
using namespace votca::csg;

class CsgStatApp : public CsgApplication {
 public:
  CsgStatApp() : _write_every(0) {}

  string ProgramName() override { return "csg_partial_rdf"; }
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

  CsgApplication::Worker *ForkWorker() override {
    return _rdf_calculator.ForkWorker();
  }

  void MergeWorker(CsgApplication::Worker *worker) override {
    _rdf_calculator.MergeWorker(worker);
  }

 public:
  RDFCalculator _rdf_calculator;
  Index _write_every;
};

void CsgStatApp::HelpText(ostream &out) {
  out << "Calculate spatially confined rdfs";
}

void CsgStatApp::Initialize() {
  CsgApplication::Initialize();
  AddProgramOptions("Specific options")("options",
                                        boost::program_options::value<string>(),
                                        "  options file defining the rdfs")(
      "subvolume_radius", boost::program_options::value<double>(),
      "Rdf calc. in spherical subvolume of this radius (from center of box)")(
      "do-vol-corr", "Correct for subvolume truncation in rdf")(
      "write-every", boost::program_options::value<Index>(&_write_every),
      " (UNIMPLEMENTED) write after every block of this length, "
      "if --blocking   is set, the averages are cleared after every output")(
      "do-blocks", "  write output for blocking analysis");
}

bool CsgStatApp::EvaluateOptions() {
  CsgApplication::EvaluateOptions();
  CheckRequired("options");
  CheckRequired("subvolume_radius");
  CheckRequired("trj", "no trajectory file specified");

  _rdf_calculator.LoadOptions(OptionsMap()["options"].as<string>());

  _rdf_calculator.SetSubvolRadius(
      OptionsMap()["subvolume_radius"].as<double>());

  _rdf_calculator.WriteEvery(_write_every);
  if (OptionsMap().count("do-blocks")) {
    _rdf_calculator.DoBlocks(true);
  }

  if (OptionsMap().count("do-vol-corr")) {
    _rdf_calculator.DoVolumeCorrection(true);
  }

  _rdf_calculator.Initialize();
  return true;
}

void CsgStatApp::BeginEvaluate(Topology *top, Topology *top_ref) {
  _rdf_calculator.BeginEvaluate(top, top_ref);
}

void CsgStatApp::EndEvaluate() { _rdf_calculator.EndEvaluate(); }

int main(Index argc, char **argv) {
  CsgStatApp app;
  app.Exec(argc, argv);
}
