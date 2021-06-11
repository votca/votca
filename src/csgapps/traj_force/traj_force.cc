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
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

// VOTCA includes
#include <votca/tools/linalg.h>
#include <votca/tools/table.h>

// Local VOTCA includes
#include <votca/csg/beadlist.h>

// Local private includes
#include "traj_force.h"

int main(int argc, char **argv) {
  TrajForce app;
  return app.Exec(argc, argv);
}

void TrajForce::Initialize(void) {
  CsgApplication::Initialize();
  AddProgramOptions()(
      "scale",
      boost::program_options::value<double>(&scale_)->default_value(-1.0),
      "  scaling factor for trajectory forces")(
      "trj-force", boost::program_options::value<string>(),
      "  atomistic reference "
      "trajectory containing forces "
      "to add/subtract")("out", boost::program_options::value<string>(),
                         "  output "
                         "trajectory "
                         "file with "
                         "resultant "
                         "forces");
}

bool TrajForce::EvaluateOptions() {
  CsgApplication::EvaluateOptions();
  CheckRequired("trj", "no trajectory file specified");
  CheckRequired("trj-force", "no reference trajectory file specified");
  CheckRequired("out", "no output trajectory file specified");

  return true;
}

void TrajForce::BeginEvaluate(Topology *top, Topology *) {
  top_force_.CopyTopologyData(top);
  trjreader_force_ =
      TrjReaderFactory().Create(OptionsMap()["trj-force"].as<string>());
  if (trjreader_force_ == nullptr) {
    throw runtime_error(string("input format not supported: ") +
                        OptionsMap()["trj-force"].as<string>());
  }
  // open the trajectory
  trjreader_force_->Open(OptionsMap()["trj-force"].as<string>());
  // read in first frame
  trjreader_force_->FirstFrame(top_force_);

  // output trajectory file
  trjwriter_ = TrjWriterFactory().Create(OptionsMap()["out"].as<string>());
  if (trjwriter_ == nullptr) {
    throw runtime_error(string("output trajectory format not supported: ") +
                        OptionsMap()["out"].as<string>());
  }
  bool append = true;
  trjwriter_->Open(OptionsMap()["out"].as<string>(), append);
}

void TrajForce::EndEvaluate() {
  cout << "\nWe are done, thank you very much!" << endl;
  trjreader_force_->Close();
  trjwriter_->Close();
}

void TrajForce::WriteOutFiles() {}

void TrajForce::EvalConfiguration(Topology *conf, Topology *) {
  if (conf->BeadCount() != top_force_.BeadCount()) {
    throw std::runtime_error(
        "number of beads in topology and reference force topology does not "
        "match");
  }
  for (votca::Index i = 0; i < conf->BeadCount(); ++i) {

    // \todo check why "conf" HasForce() is false
    // Since "conf" topology Force is set to false
    // for now using  top_force_ to store resultant output forces

    top_force_.getBead(i)->F() =
        conf->getBead(i)->getF() + scale_ * top_force_.getBead(i)->getF();
    Eigen::Vector3d d =
        conf->getBead(i)->getPos() - top_force_.getBead(i)->getPos();
    if (d.norm() > 1e-6) {
      throw std::runtime_error(
          "One or more bead positions in trajectory and reference force "
          "trajectory differ by more than 1e-6");
    }
  }

  trjwriter_->Write(&top_force_);
  trjreader_force_->NextFrame(top_force_);
}
