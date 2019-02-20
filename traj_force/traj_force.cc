/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include "traj_force.h"
#include <boost/numeric/ublas/vector.hpp>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <votca/csg/beadlist.h>
#include <votca/tools/linalg.h>
#include <votca/tools/table.h>

int main(int argc, char **argv) {
  TrajForce app;
  return app.Exec(argc, argv);
}

void TrajForce::Initialize(void) {
  CsgApplication::Initialize();
  AddProgramOptions()(
      "scale",
      boost::program_options::value<double>(&_scale)->default_value(-1.0),
      "  scaling factor for trajectory forces")(
      "trj-force", boost::program_options::value<string>(),
      "  atomistic reference trajectory containing forces to add/subtract")(
      "out", boost::program_options::value<string>(),
      "  output trajectory file with resultant forces");
}

bool TrajForce::EvaluateOptions() {
  CsgApplication::EvaluateOptions();
  CheckRequired("trj", "no trajectory file specified");
  CheckRequired("trj-force", "no reference trajectory file specified");
  CheckRequired("out", "no output trajectory file specified");

  return true;
}

void TrajForce::BeginEvaluate(Topology *top, Topology *top_atom) {
  _top_force.CopyTopologyData(top);
  _trjreader_force =
      TrjReaderFactory().Create(_op_vm["trj-force"].as<string>());
  if (_trjreader_force == NULL)
    throw runtime_error(string("input format not supported: ") +
                        _op_vm["trj-force"].as<string>());
  // open the trajectory
  _trjreader_force->Open(_op_vm["trj-force"].as<string>());
  // read in first frame
  _trjreader_force->FirstFrame(_top_force);

  // output trajectory file
  _trjwriter = TrjWriterFactory().Create(_op_vm["out"].as<string>());
  if (_trjwriter == NULL)
    throw runtime_error(string("output trajectory format not supported: ") +
                        _op_vm["out"].as<string>());
  bool append = true;
  _trjwriter->Open(_op_vm["out"].as<string>(), append);
}

void TrajForce::EndEvaluate() {
  cout << "\nWe are done, thank you very much!" << endl;
  _trjreader_force->Close();
  _trjwriter->Close();
  delete _trjreader_force;
  delete _trjwriter;
}

void TrajForce::WriteOutFiles() {}

void TrajForce::EvalConfiguration(Topology *conf, Topology *conf_atom) {
  if (conf->BeadCount() != _top_force.BeadCount())
    throw std::runtime_error(
        "number of beads in topology and reference force topology does not "
        "match");
  for (int i = 0; i < conf->BeadCount(); ++i) {
    // conf->getBead(i)->F() += _scale*_top_force.getBead(i)->getF();

    // \todo check why "conf" HasForce() is false
    // Since "conf" topology Force is set to false
    // for now using _top_force to store resultant output forces

    _top_force.getBead(i)->F() =
        conf->getBead(i)->getF() + _scale * _top_force.getBead(i)->getF();
    Eigen::Vector3d d = conf->getBead(i)->getPos() - _top_force.getBead(i)->getPos();
    if (abs(d) > 1e-6)
      throw std::runtime_error(
          "One or more bead positions in trajectory and reference force "
          "trajectory differ by more than 1e-6");
  }
  //  cout << conf->HasForce() << endl;
  //  cout << _top_force.HasForce() << endl;
  _trjwriter->Write(&_top_force);
  _trjreader_force->NextFrame(_top_force);
}
