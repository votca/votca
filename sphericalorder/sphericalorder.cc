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

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>
#include <math.h>

#include <stdlib.h>
#include <votca/csg/cgengine.h>
#include <votca/csg/csgapplication.h>
#include <votca/tools/average.h>
#include <votca/tools/tokenizer.h>

using namespace std;
using namespace votca::csg;

class CGOrderParam : public CsgApplication {
 public:
  string ProgramName() override { return "sphericalorder"; }

  void HelpText(ostream &out) override {

    out << "!! EXPERIMENTAL !! Calculate spherical order parameter.\n"
           " Needs non-spherical beads in mapping.\n\n";
  }

  void Initialize() override {
    CsgApplication::Initialize();
    AddProgramOptions()(
        "filter",
        boost::program_options::value<string>(&_filter)->default_value("*"),
        "filter molecule names")(
        "radialcut",
        boost::program_options::value<double>(&_radialcutoff)
            ->default_value(0.0),
        "radial cutoff: distance from center where bead is considered")(
        "minrad",
        boost::program_options::value<double>(&_minrad)->default_value(0.0),
        "minimal distance a parcle has to be apart from center to be "
        "considerd")(
        "refmol",
        boost::program_options::value<string>(&_refmol)->default_value(""),
        "Reference molecule")(
        "rbinw",
        boost::program_options::value<double>(&_rbinw)->default_value(0),
        "Do multiple r_bins multiple histograms");
  }

  bool EvaluateOptions() override {
    CsgApplication::EvaluateOptions();
    // CheckRequired("radialcut");
    return true;
  }

  bool DoTrajectory() override { return true; }
  bool DoMapping() override { return true; }

  void BeginEvaluate(Topology *top, Topology *) override {

    string filter;

    filter = OptionsMap()["filter"].as<string>();

    _minrad = 0;

    _radialcutoff = OptionsMap()["radialcut"].as<double>();
    _minrad = OptionsMap()["minrad"].as<double>();
    _refmol = OptionsMap()["refmol"].as<string>();
    _rbinw = OptionsMap()["rbinw"].as<double>();

    if (_rbinw == 0 && _radialcutoff <= 0) {
      throw runtime_error("_radialcut > 0 has to be specified");
    }

    setFilter(filter);

    _file_u.open("hist_u.xvg");
    if (!_file_u) {
      throw runtime_error("cannot open hist_u.xvg for output");
    }
    _file_v.open("hist_v.xvg");
    if (!_file_v) {
      throw runtime_error("cannot open hist_v.xvg for output");
    }
    _file_w.open("hist_w.xvg");
    if (!_file_w) {
      throw runtime_error("cannot open hist_w.xvg for output");
    }

    _n = 0;

    Eigen::Matrix3d box = top->getBox();
    Eigen::Vector3d a = box.col(0);

    if (_refmol == "") {

      _ref = box.rowwise().sum() / 2;

      cout << "Refernce is center of box " << _ref << endl;
    }

    boxl = a.norm() / 2;
    if (_rbinw > 0) {
      _rbins = (votca::Index)(boxl / _rbinw) + 1;
      cout << "radial bins " << _rbins << endl;
    } else {
      _rbins = 1;
      cout << "considering atoms between " << _minrad << " and "
           << _radialcutoff << endl;
    }

    _nbin = 100;
    _hist_u = Eigen::MatrixXd::Zero(_rbins, _nbin);
    _hist_v = Eigen::MatrixXd::Zero(_rbins, _nbin);
    _hist_w = Eigen::MatrixXd::Zero(_rbins, _nbin);
    _nmol = Eigen::VectorXi::Zero(_rbins);
  }

  void EndEvaluate() override {

    cout << "Average number of molecules within cutoff " << endl;
    for (votca::Index i = 0; i < _rbins; i++) {
      cout << (double)i * _rbinw << " " << (double)_nmol[i] / (double)_n << endl;
    }

    double exp_value = 1.0 / (double)_nbin;
    double orderparam = 0;

    for (votca::Index n = 0; n < _nbin; n++) {
      _hist_u(0, n) /= (double)_nmol[0];  // normalize to numberframes and avg.
                                          // number of molecules
      _hist_v(0, n) /= (double)_nmol[0];
      _hist_w(0, n) /= (double)_nmol[0];

      _file_u << (double)n * 2 / double(_nbin-1) << " " << _hist_u(0, n) << endl;
      _file_v << (double)n * 2 / double(_nbin-1) << " " << _hist_v(0, n) << endl;
      _file_w << (double)n * 2 / double(_nbin-1) << " " << _hist_w(0, n) << endl;

      orderparam += (_hist_u(0, n) - exp_value) * (_hist_u(0, n) - exp_value);
    }

    orderparam = sqrt(orderparam / (double)_nbin);

    cout << "Orderparam " << _radialcutoff << " " << orderparam << endl;

    _file_u.close();
    _file_v.close();
    _file_w.close();
  };

  void EvalConfiguration(Topology *conf, Topology * = nullptr) override {

    Eigen::Vector3d eR;
    votca::Index nu, nv, nw;
    Eigen::Vector3d u, v, w;

    if (_refmol != "") {
      for (Bead *bead : conf->Beads()) {
        if (votca::tools::wildcmp(_refmol.c_str(), bead->getName().c_str())) {
          _ref = bead->getPos();
        }
      }
    }

    for (Bead *bead : conf->Beads()) {
      if (!votca::tools::wildcmp(_filter.c_str(), bead->getName().c_str())) {
        continue;
      }
      if (votca::tools::wildcmp(_refmol.c_str(), bead->getName().c_str())) {
        continue;
      }

      eR = bead->getPos() - _ref;
      if ((eR.norm() < _radialcutoff && eR.norm() > _minrad) || _rbins != 1) {
        // cout << eR << endl;
        votca::Index rb = 0;
        if (_rbinw > 0) {
          rb = (votca::Index)((eR.norm()) / boxl * (double)_rbins);
        }
        if (rb >= _rbins) {
          continue;
        }

        eR.normalize();
        u = bead->getU();
        v = bead->getV();
        w = bead->getW();
        u.normalize();
        v.normalize();
        w.normalize();

        nu = (votca::Index)((eR.dot(u) + 1) / 2) * _nbin;
        nv = (votca::Index)((eR.dot(v) + 1) / 2) * _nbin;
        nw = (votca::Index)((eR.dot(w) + 1) / 2) * _nbin;

        _hist_u(rb, nu) += 1;
        _hist_v(rb, nv) += 1;
        _hist_w(rb, nw) += 1;
        _nmol[rb]++;
      }
    }

    _n++;
  }

  void setOut(const string &filename) { _filename = filename; }

  void setFilter(const string &filter) { _filter = filter; }

 protected:
  ofstream _file;
  string _filename;
  votca::Index _n;
  Eigen::Vector3d _ref;
  ofstream _file_u;
  ofstream _file_v;
  ofstream _file_w;
  Eigen::MatrixXd _hist_u;
  Eigen::MatrixXd _hist_v;
  Eigen::MatrixXd _hist_w;
  votca::Index _nbin;
  Eigen::VectorXi _nmol;
  double _radialcutoff;
  double _minrad;
  votca::Index _rbins;
  double _rbinw;
  double boxl;

  string _filter;
  string _refmol;
};

int main(int argc, char **argv) {
  CGOrderParam app;
  return app.Exec(argc, argv);
}
