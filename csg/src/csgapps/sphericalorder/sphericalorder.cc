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
#include <cstdlib>
#include <fstream>
#include <iostream>

// Third party includes
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

// VOTCA includes
#include <votca/tools/average.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include <votca/csg/cgengine.h>
#include <votca/csg/csgapplication.h>

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
        boost::program_options::value<string>(&filter_)->default_value("*"),
        "filter molecule names")(
        "radialcut",
        boost::program_options::value<double>(&radialcutoff_)
            ->default_value(0.0),
        "radial cutoff: distance from center where bead is considered")(
        "minrad",
        boost::program_options::value<double>(&minrad_)->default_value(0.0),
        "minimal distance a parcle has to be apart from center to be "
        "considerd")(
        "refmol",
        boost::program_options::value<string>(&refmol_)->default_value(""),
        "Reference molecule")(
        "rbinw",
        boost::program_options::value<double>(&rbinw_)->default_value(0),
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

    minrad_ = 0;

    radialcutoff_ = OptionsMap()["radialcut"].as<double>();
    minrad_ = OptionsMap()["minrad"].as<double>();
    refmol_ = OptionsMap()["refmol"].as<string>();
    rbinw_ = OptionsMap()["rbinw"].as<double>();

    if (rbinw_ == 0 && radialcutoff_ <= 0) {
      throw runtime_error(" radialcut_ > 0 has to be specified");
    }

    setFilter(filter);

    file_u_.open("hist_u.xvg");
    if (!file_u_) {
      throw runtime_error("cannot open hist_u.xvg for output");
    }
    file_v_.open("hist_v.xvg");
    if (!file_v_) {
      throw runtime_error("cannot open hist_v.xvg for output");
    }
    file_w_.open("hist_w.xvg");
    if (!file_w_) {
      throw runtime_error("cannot open hist_w.xvg for output");
    }

    n_ = 0;

    Eigen::Matrix3d box = top->getBox();
    Eigen::Vector3d a = box.col(0);

    if (refmol_ == "") {

      ref_ = box.rowwise().sum() / 2;

      cout << "Refernce is center of box " << ref_ << endl;
    }

    boxl = a.norm() / 2;
    if (rbinw_ > 0) {
      rbins_ = (votca::Index)(boxl / rbinw_) + 1;
      cout << "radial bins " << rbins_ << endl;
    } else {
      rbins_ = 1;
      cout << "considering atoms between " << minrad_ << " and "
           << radialcutoff_ << endl;
    }

    nbin_ = 100;
    hist_u_ = Eigen::MatrixXd::Zero(rbins_, nbin_);
    hist_v_ = Eigen::MatrixXd::Zero(rbins_, nbin_);
    hist_w_ = Eigen::MatrixXd::Zero(rbins_, nbin_);
    nmol_ = Eigen::VectorXi::Zero(rbins_);
  }

  void EndEvaluate() override {

    cout << "Average number of molecules within cutoff " << endl;
    for (votca::Index i = 0; i < rbins_; i++) {
      cout << (double)i * rbinw_ << " " << (double)nmol_[i] / (double)n_
           << endl;
    }

    double exp_value = 1.0 / (double)nbin_;
    double orderparam = 0;

    for (votca::Index n = 0; n < nbin_; n++) {
      hist_u_(0, n) /= (double)nmol_[0];  // normalize to numberframes and avg.
                                          // number of molecules
      hist_v_(0, n) /= (double)nmol_[0];
      hist_w_(0, n) /= (double)nmol_[0];

      file_u_ << (double)n * 2 / double(nbin_ - 1) << " " << hist_u_(0, n)
              << endl;
      file_v_ << (double)n * 2 / double(nbin_ - 1) << " " << hist_v_(0, n)
              << endl;
      file_w_ << (double)n * 2 / double(nbin_ - 1) << " " << hist_w_(0, n)
              << endl;

      orderparam += (hist_u_(0, n) - exp_value) * (hist_u_(0, n) - exp_value);
    }

    orderparam = sqrt(orderparam / (double)nbin_);

    cout << "Orderparam " << radialcutoff_ << " " << orderparam << endl;

    file_u_.close();
    file_v_.close();
    file_w_.close();
  };

  void EvalConfiguration(Topology *conf, Topology * = nullptr) override {

    Eigen::Vector3d eR;
    votca::Index nu, nv, nw;
    Eigen::Vector3d u, v, w;

    if (refmol_ != "") {
      for (const auto &bead : conf->Beads()) {
        if (votca::tools::wildcmp(refmol_, bead.getName())) {
          ref_ = bead.getPos();
        }
      }
    }

    for (const auto &bead : conf->Beads()) {
      if (!votca::tools::wildcmp(filter_, bead.getName())) {
        continue;
      }
      if (votca::tools::wildcmp(refmol_, bead.getName())) {
        continue;
      }

      eR = bead.getPos() - ref_;
      if ((eR.norm() < radialcutoff_ && eR.norm() > minrad_) || rbins_ != 1) {
        // cout << eR << endl;
        votca::Index rb = 0;
        if (rbinw_ > 0) {
          rb = (votca::Index)((eR.norm()) / boxl * (double)rbins_);
        }
        if (rb >= rbins_) {
          continue;
        }

        eR.normalize();
        u = bead.getU();
        v = bead.getV();
        w = bead.getW();
        u.normalize();
        v.normalize();
        w.normalize();

        nu = (votca::Index)((eR.dot(u) + 1) / 2) * nbin_;
        nv = (votca::Index)((eR.dot(v) + 1) / 2) * nbin_;
        nw = (votca::Index)((eR.dot(w) + 1) / 2) * nbin_;

        hist_u_(rb, nu) += 1;
        hist_v_(rb, nv) += 1;
        hist_w_(rb, nw) += 1;
        nmol_[rb]++;
      }
    }

    n_++;
  }

  void setOut(const string &filename) { filename_ = filename; }

  void setFilter(const string &filter) { filter_ = filter; }

 protected:
  ofstream file_;
  string filename_;
  votca::Index n_;
  Eigen::Vector3d ref_;
  ofstream file_u_;
  ofstream file_v_;
  ofstream file_w_;
  Eigen::MatrixXd hist_u_;
  Eigen::MatrixXd hist_v_;
  Eigen::MatrixXd hist_w_;
  votca::Index nbin_;
  Eigen::VectorXi nmol_;
  double radialcutoff_;
  double minrad_;
  votca::Index rbins_;
  double rbinw_;
  double boxl;

  string filter_;
  string refmol_;
};

int main(int argc, char **argv) {
  CGOrderParam app;
  return app.Exec(argc, argv);
}
