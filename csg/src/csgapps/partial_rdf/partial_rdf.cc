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

// Standard includes
#include <fstream>
#include <iostream>

// Third party includes
#include <boost/program_options.hpp>

// Local VOTCA includes
#include <votca/csg/cgengine.h>
#include <votca/csg/version.h>

// Local internal includes
#include "rdf_calculator.h"

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

#include <cstdlib>
#include <votca/csg/csgapplication.h>

// using namespace votca::tools;
using namespace std;
using namespace votca::csg;

class CsgPartialRdfApp : public CsgApplication {
 public:
  CsgPartialRdfApp() : write_every_(0) {}

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

  std::unique_ptr<CsgApplication::Worker> ForkWorker() override {
    return rdf_calculator_.ForkWorker();
  }

  void MergeWorker(CsgApplication::Worker *worker) override {
    rdf_calculator_.MergeWorker(worker);
  }

 public:
  RDFCalculator rdf_calculator_;
  votca::Index write_every_;
};

void CsgPartialRdfApp::HelpText(ostream &out) {
  out << "Calculate spatially confined rdfs";
}

void CsgPartialRdfApp::Initialize() {
  CsgApplication::Initialize();
  AddProgramOptions("Specific options")("options",
                                        boost::program_options::value<string>(),
                                        "  options file defining the rdfs")(
      "subvolume_radius", boost::program_options::value<double>(),
      "Rdf calc. in spherical subvolume of this radius (from center of box)")(
      "do-vol-corr", "Correct for subvolume truncation in rdf")(
      "write-every", boost::program_options::value<votca::Index>(&write_every_),
      " (UNIMPLEMENTED) write after every block of this length, "
      "if --blocking   is set, the averages are cleared after every output")(
      "do-blocks", "  write output for blocking analysis");
}

bool CsgPartialRdfApp::EvaluateOptions() {
  CsgApplication::EvaluateOptions();
  CheckRequired("options");
  CheckRequired("subvolume_radius");
  CheckRequired("trj", "no trajectory file specified");

  rdf_calculator_.LoadOptions(OptionsMap()["options"].as<string>());

  rdf_calculator_.SetSubvolRadius(
      OptionsMap()["subvolume_radius"].as<double>());

  rdf_calculator_.WriteEvery(write_every_);
  if (OptionsMap().count("do-blocks")) {
    rdf_calculator_.DoBlocks(true);
  }

  if (OptionsMap().count("do-vol-corr")) {
    rdf_calculator_.DoVolumeCorrection(true);
  }

  rdf_calculator_.Initialize();
  return true;
}

void CsgPartialRdfApp::BeginEvaluate(Topology *top, Topology *top_ref) {
  rdf_calculator_.BeginEvaluate(top, top_ref);
}

void CsgPartialRdfApp::EndEvaluate() { rdf_calculator_.EndEvaluate(); }

int main(int argc, char **argv) {
  CsgPartialRdfApp app;
  app.Exec(argc, argv);
}
