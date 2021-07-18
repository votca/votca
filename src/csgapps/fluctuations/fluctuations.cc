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

// Third party includes
#include <boost/program_options.hpp>

// VOTCA includes
#include <votca/tools/histogramnew.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include <votca/csg/cgengine.h>
#include <votca/csg/csgapplication.h>

using namespace std;
using namespace votca::csg;

class CsgFluctuations : public CsgApplication {
  string ProgramName() override { return "fluctuations"; }
  void HelpText(ostream &out) override {
    out << "calculate density fluctuations in subvolumes of the simulation "
           "box.";
    out << "Subolumes can be either cubic slabs in dimensions (x|y|z) or "
           "spherical";
    out << "slabs with respect to either the center of box or a reference "
           "molecule";
  }

  // some program options are added here

  void Initialize() override {
    CsgApplication::Initialize();
    // add program option to pick molecule
    AddProgramOptions("Fluctuation options")(
        "filter",
        boost::program_options::value<string>(&filter_)->default_value("*"),
        "filter molecule names")("rmax",
                                 boost::program_options::value<double>(),
                                 "maximal distance to be considered")(
        "rmin",
        boost::program_options::value<double>(&rmin_)->default_value(0.0),
        "minimal distance to be considered")(
        "refmol",
        boost::program_options::value<string>(&refmol_)->default_value(""),
        "Reference molecule")(
        "nbin",
        boost::program_options::value<votca::Index>(&nbins_)->default_value(
            100),
        "Number of bins")(
        "geometry", boost::program_options::value<string>(),
        "(sphere|x|y|z) Take radial or x, y, z slabs from rmin to rmax")(
        "outfile",
        boost::program_options::value<string>(&outfilename_)
            ->default_value("fluctuations.dat"),
        "Output file");
  }
  bool EvaluateOptions() override {
    CsgApplication::EvaluateOptions();
    CheckRequired("rmax");
    CheckRequired("geometry");
    return true;
  }

  // we want to process a trajectory
  bool DoTrajectory() override { return true; }
  bool DoMapping() override { return true; }

  void BeginEvaluate(Topology *top, Topology *) override {
    filter_ = OptionsMap()["filter"].as<string>();
    refmol_ = OptionsMap()["refmol"].as<string>();
    rmin_ = OptionsMap()["rmin"].as<double>();
    rmax_ = OptionsMap()["rmax"].as<double>();
    nbins_ = OptionsMap()["nbin"].as<votca::Index>();
    outfilename_ = OptionsMap()["outfile"].as<string>();
    geometryinput_ = OptionsMap()["geometry"].as<string>();
    nframes_ = 0;

    do_spherical_ = false;

    if (geometryinput_ == "sphere") {
      cout << "Doing spherical slabs" << endl;
      do_spherical_ = true;
    } else if (geometryinput_ == "x") {
      cout << "Doing slabs along x-axis" << endl;
      dim_ = 0;
    } else if (geometryinput_ == "y") {
      cout << "Doing slabs along  y-axis" << endl;
      dim_ = 1;
    } else if (geometryinput_ == "z") {
      cout << "Doing slabs along  z-axis" << endl;
      dim_ = 2;
    } else {
      throw std::runtime_error("Unrecognized geometry option. (sphere|x|y|z)");
    }

    N_avg_ = Eigen::VectorXd::Zero(nbins_);
    N_sq_avg_ = Eigen::VectorXd::Zero(nbins_);

    if (do_spherical_) {
      cout << "Calculating fluctions for " << rmin_ << "<r<" << rmax_;
      cout << "using " << nbins_ << " bins" << endl;
    } else {
      cout << "Calculating fluctions for " << rmin_ << "<" << geometryinput_
           << "<" << rmax_;
      cout << "using " << nbins_ << " bins" << endl;
    }

    if (refmol_ == "" && do_spherical_) {
      Eigen::Matrix3d box = top->getBox();
      ref_ = box.rowwise().sum() / 2;

      cout << "Reference is center of box " << ref_ << endl;
    }

    outfile_.open(outfilename_);
    if (!outfile_) {
      throw runtime_error("cannot open" + outfilename_ + " for output");
    }
  }

  // write out results in EndEvaluate
  void EndEvaluate() override;
  // do calculation in this function
  void EvalConfiguration(Topology *conf, Topology *top_ref) override;

 protected:
  // number of particles in dV
  votca::Index nbins_;
  Eigen::VectorXd N_avg_;
  // sqare
  Eigen::VectorXd N_sq_avg_;
  string filter_;
  string refmol_;
  double rmax_;
  double rmin_;
  Eigen::Vector3d ref_;
  votca::Index nframes_;
  string outfilename_;
  ofstream outfile_;
  string geometryinput_;
  bool do_spherical_;
  votca::Index dim_;
};

int main(int argc, char **argv) {
  CsgFluctuations app;

  return app.Exec(argc, argv);
}

void CsgFluctuations::EvalConfiguration(Topology *conf, Topology *) {

  if (refmol_ != "") {
    for (const auto &bead : conf->Beads()) {
      if (votca::tools::wildcmp(refmol_, bead.getName())) {
        ref_ = bead.getPos();
        cout << " Solute pos " << ref_ << endl;
      }
    }
  }

  votca::tools::HistogramNew hist;
  hist.Initialize(rmin_, rmax_, nbins_);

  /* check how many molecules are in each bin*/
  for (const auto &bead : conf->Beads()) {
    if (!votca::tools::wildcmp(filter_, bead.getName())) {
      continue;
    }
    double r = 0;
    if (do_spherical_) {
      Eigen::Vector3d eR = bead.getPos() - ref_;
      r = eR.norm();
    } else {
      Eigen::Vector3d eR = bead.getPos();
      if (dim_ == 0) {
        r = eR.x();
      } else if (dim_ == 1) {
        r = eR.y();
      } else if (dim_ == 2) {
        r = eR.z();
      }
    }
    hist.Process(r);
  }

  /* update averages*/
  N_avg_ += hist.data().y();
  N_sq_avg_ += hist.data().y().cwiseAbs2();

  nframes_++;
}

// output everything when processing frames is done
void CsgFluctuations::EndEvaluate() {
  cout << "Writing results to " << outfilename_ << endl;
  outfile_ << "# radius number_fluct avg_number" << endl;

  for (votca::Index i = 0; i < nbins_; i++) {
    N_avg_[i] /= (double)nframes_;
    N_sq_avg_[i] /= (double)nframes_;
  }
  for (votca::Index i = 0; i < nbins_; i++) {
    outfile_ << rmin_ + (double)i * (rmax_ - rmin_) / (double)nbins_ << " ";
    outfile_ << (N_sq_avg_[i] - N_avg_[i] * N_avg_[i]) / N_avg_[i]
             << " ";  // fluctuation
    outfile_ << N_avg_[i] << endl;
  }
}

// add our user program options
