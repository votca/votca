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

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <votca/csg/cgengine.h>
#include <votca/csg/csgapplication.h>
#include <votca/tools/average.h>
#include <votca/tools/tokenizer.h>

// using namespace votca::tools;
using namespace std;
using namespace votca::csg;

class CsgFluctuations : public CsgApplication {
  string ProgramName() { return "fluctuations"; }
  void HelpText(ostream &out) {
    out << "calculate density fluctuations in subvolumes of the simulation "
           "box.";
    out << "Subolumes can be either cubic slabs in dimensions (x|y|z) or "
           "spherical";
    out << "slabs with respect to either the center of box or a reference "
           "molecule";
  }

  // some program options are added here

  void Initialize() {
    CsgApplication::Initialize();
    // add program option to pick molecule
    AddProgramOptions("Fluctuation options")(
        "filter",
        boost::program_options::value<string>(&_filter)->default_value("*"),
        "filter molecule names")("rmax",
                                 boost::program_options::value<double>(),
                                 "maximal distance to be considered")(
        "rmin",
        boost::program_options::value<double>(&_rmin)->default_value(0.0),
        "minimal distance to be considered")(
        "refmol",
        boost::program_options::value<string>(&_refmol)->default_value(""),
        "Reference molecule")(
        "nbin", boost::program_options::value<int>(&_nbins)->default_value(100),
        "Number of bins")(
        "geometry", boost::program_options::value<string>(),
        "(sphere|x|y|z) Take radial or x, y, z slabs from rmin to rmax")(
        "outfile",
        boost::program_options::value<string>(&_outfilename)
            ->default_value("fluctuations.dat"),
        "Output file");
  }
  bool EvaluateOptions() {
    CsgApplication::EvaluateOptions();
    CheckRequired("rmax");
    CheckRequired("geometry");
    return true;
  }

  // we want to process a trajectory
  bool DoTrajectory() { return true; }
  bool DoMapping() { return true; }

  void BeginEvaluate(Topology *top, Topology *top_atom) {
    _filter = OptionsMap()["filter"].as<string>();
    _refmol = OptionsMap()["refmol"].as<string>();
    _rmin = OptionsMap()["rmin"].as<double>();
    _rmax = OptionsMap()["rmax"].as<double>();
    _nbins = OptionsMap()["nbin"].as<int>();
    _outfilename = OptionsMap()["outfile"].as<string>();
    _geometryinput = OptionsMap()["geometry"].as<string>();
    _nframes = 0;

    _do_spherical = false;

    if (_geometryinput == "sphere") {
      cout << "Doing spherical slabs" << endl;
      _do_spherical = true;
    } else if (_geometryinput == "x") {
      cout << "Doing slabs along x-axis" << endl;
      _dim = 0;
    } else if (_geometryinput == "y") {
      cout << "Doing slabs along  y-axis" << endl;
      _dim = 1;
    } else if (_geometryinput == "z") {
      cout << "Doing slabs along  z-axis" << endl;
      _dim = 2;
    } else {
      cout << "Unrecognized geometry option. (sphere|x|y|z)" << endl;
      exit(0);
    }

    _N_avg = new double[_nbins];
    _N_sq_avg = new double[_nbins];
    N = new int[_nbins];
    for (int i = 0; i < _nbins; i++) {
      _N_avg[i] = 0;
      _N_sq_avg[i] = 0;
    }

    if (_do_spherical) {
      cout << "Calculating fluctions for " << _rmin << "<r<" << _rmax;
      cout << "using " << _nbins << " bins" << endl;
    } else {
      cout << "Calculating fluctions for " << _rmin << "<" << _geometryinput
           << "<" << _rmax;
      cout << "using " << _nbins << " bins" << endl;
    }

    if (_refmol == "" && _do_spherical) {
      Eigen::Matrix3d box = top->getBox();
      Eigen::Vector3d a = box.col(0);
      Eigen::Vector3d b = box.col(1);
      Eigen::Vector3d c = box.col(2);
      _ref = (a + b + c) / 2;

      cout << "Refernce is center of box " << _ref << endl;
    }

    _outfile.open(_outfilename.c_str());
    if (!_outfile) throw runtime_error("cannot open outfile for output");
  }

  // write out results in EndEvaluate
  void EndEvaluate();
  // do calculation in this function
  void EvalConfiguration(Topology *top, Topology *top_ref);

 protected:
  // number of particles in dV
  int _nbins;
  double *_N_avg;
  // sqare
  double *_N_sq_avg;
  int *N;
  string _filter;
  string _refmol;
  double _rmax;
  double _rmin;
  Eigen::Vector3d _ref;
  int _nframes;
  string _outfilename;
  ofstream _outfile;
  string _geometryinput;
  bool _do_spherical;
  int _dim;
};

int main(int argc, char **argv) {
  CsgFluctuations app;

  return app.Exec(argc, argv);
}

void CsgFluctuations::EvalConfiguration(Topology *conf,
                                        Topology *conf_atom = 0) {
  Eigen::Vector3d eR;
  double r = 0;
  int rbin;

  if (_refmol != "") {
    for (Bead *bead : conf->Beads()) {
      if (votca::tools::wildcmp(_refmol.c_str(), bead->getName().c_str())) {
        _ref = bead->getPos();
        cout << " Solute pos " << _ref << endl;
      }
    }
  }

  for (int i = 0; i < _nbins; i++) {
    N[i] = 0;
  }

  /* check how many molecules are in each bin*/
  for (Bead *bead : conf->Beads()) {
    if (!votca::tools::wildcmp(_filter.c_str(), bead->getName().c_str()))
      continue;

    if (_do_spherical) {
      eR = bead->getPos() - _ref;
      r = eR.norm();
    } else {
      eR = bead->getPos();
      if (_dim == 0)
        r = eR.x();
      else if (_dim == 1)
        r = eR.y();
      else if (_dim == 2)
        r = eR.z();
    }
    if (r > _rmin && r < _rmax) {
      rbin = (int)_nbins * (double)((r - _rmin) / (_rmax - _rmin));
      N[rbin]++;
    }
  }

  /* update averages*/
  for (int i = 0; i < _nbins; i++) {
    _N_avg[i] += N[i];
    _N_sq_avg[i] += N[i] * N[i];
  }

  _nframes++;
}

// output everything when processing frames is done
void CsgFluctuations::EndEvaluate() {
  cout << "Writing results to " << _outfilename << endl;
  _outfile << "# radius number_fluct avg_number" << endl;

  for (int i = 0; i < _nbins; i++) {
    _N_avg[i] /= _nframes;
    _N_sq_avg[i] /= _nframes;
  }
  for (int i = 0; i < _nbins; i++) {
    _outfile << _rmin + i * (_rmax - _rmin) / _nbins << " ";
    _outfile << (_N_sq_avg[i] - _N_avg[i] * _N_avg[i]) / _N_avg[i]
             << " ";  // fluctuation
    _outfile << _N_avg[i] << endl;
  }
}

// add our user program options
