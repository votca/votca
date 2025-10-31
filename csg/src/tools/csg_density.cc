/*
 * Copyright 2009-2025 The VOTCA Development Team (http://www.votca.org)
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

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/histogram.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/csg/csgapplication.h"

using namespace std;
using namespace votca::csg;

class CsgDensityApp : public CsgApplication {
  string ProgramName() override { return "csg_density"; }
  void HelpText(ostream &out) override {
    out << "Calculates the mass density distribution along a box axis or "
           "radial density profile from reference point";
  }

  // some program options are added here
  void Initialize() override;

  // we want to process a trajectory
  bool DoTrajectory() override { return true; }
  bool DoMapping() override { return true; }
  bool DoMappingDefault(void) override { return false; }

  // write out results in EndEvaluate
  void EndEvaluate() override;
  void BeginEvaluate(Topology *top, Topology *top_atom) override;
  void EvalConfiguration(Topology *top, Topology *top_ref) override;

  bool EvaluateOptions() override {
    CsgApplication::EvaluateOptions();
    CheckRequired("out", "no output topology specified");
    CheckRequired("trj", "no trajectory file specified");
    return true;
  };

 protected:
  string filter_, out_;
  votca::tools::Histogram dist_;
  string dens_type_;
  double rmax_;
  votca::Index nbin_;
  double scale_;
  double step_;
  votca::Index frames_;
  votca::Index nblock_;
  votca::Index block_length_;
  Eigen::Vector3d ref_;
  Eigen::Vector3d axis_;
  string axisname_;
  string molname_;
  double area_;
  void WriteDensity(votca::Index nframes, const string &suffix = "");
};

int main(int argc, char **argv) {
  CsgDensityApp app;
  return app.Exec(argc, argv);
}

void CsgDensityApp::BeginEvaluate(Topology *top, Topology *) {

  Eigen::Matrix3d box = top->getBox();
  Eigen::Vector3d a = box.col(0);
  Eigen::Vector3d b = box.col(1);
  Eigen::Vector3d c = box.col(2);

  dist_.setPeriodic(true);
  axis_ = Eigen::Vector3d::Zero();
  area_ = 0;
  if (axisname_ == "x") {
    axis_ = Eigen::Vector3d::UnitX();
    rmax_ = a.norm();
    area_ = b.cross(c).norm();
  } else if (axisname_ == "y") {
    axis_ = Eigen::Vector3d::UnitY();
    rmax_ = b.norm();
    area_ = a.cross(c).norm();
  } else if (axisname_ == "z") {
    axis_ = Eigen::Vector3d::UnitZ();
    rmax_ = c.norm();
    area_ = a.cross(b).norm();
  } else if (axisname_ == "r") {
    dist_.setPeriodic(false);
    rmax_ = min(min((a / 2).norm(), (b / 2).norm()), (c / 2).norm());
  } else {
    throw std::runtime_error("unknown axis type");
  }

  if (OptionsMap().count("rmax")) {
    rmax_ = OptionsMap()["rmax"].as<double>();
  }

  if (OptionsMap().count("block-length")) {
    block_length_ = OptionsMap()["block-length"].as<votca::Index>();
  } else {
    block_length_ = 0;
  }

  if (axisname_ == "r") {
    if (!OptionsMap().count("ref")) {
      ref_ = a / 2 + b / 2 + c / 2;
    }
    cout << "Using referece point: " << ref_ << endl;
  } else if (OptionsMap().count("ref")) {
    throw std::runtime_error(
        "reference center can only be used in case of spherical density");
  }

  nbin_ = (votca::Index)floor(rmax_ / step_);
  dist_.Initialize(0, rmax_, nbin_);

  cout << "rmax: " << rmax_ << endl;
  cout << "axis: " << axisname_ << endl;
  cout << "Bins: " << nbin_ << endl;
  frames_ = 0;
  nblock_ = 0;
}

void CsgDensityApp::EvalConfiguration(Topology *top, Topology *) {
  // loop over all molecules
  bool did_something = false;
  for (const auto &mol : top->Molecules()) {
    if (!votca::tools::wildcmp(molname_, mol.getName())) {
      continue;
    }
    votca::Index N = mol.BeadCount();
    for (votca::Index i = 0; i < N; i++) {
      const Bead *b = mol.getBead(i);
      if (!votca::tools::wildcmp(filter_, b->getName())) {
        continue;
      }
      double r;
      if (axisname_ == "r") {
        r = (top->BCShortestConnection(ref_, b->getPos()).norm());
      } else {
        r = b->getPos().dot(axis_);
      }
      if (dens_type_ == "mass") {
        dist_.Process(r, b->getMass());
      } else {
        dist_.Process(r, 1.0);
      }
      did_something = true;
    }
  }
  frames_++;
  if (!did_something) {
    throw std::runtime_error("No molecule in selection");
  }
  if (block_length_ != 0) {
    if ((nframes_ % block_length_) == 0) {
      nblock_++;
      string suffix = string("_") + boost::lexical_cast<string>(nblock_);
      WriteDensity(block_length_, suffix);
      dist_.Clear();
    }
  }
}

// output everything when processing frames is done
void CsgDensityApp::WriteDensity(votca::Index nframes, const string &suffix) {
  if (axisname_ == "r") {
    dist_.data().y() =
        scale_ /
        (double(nframes) * rmax_ / (double)nbin_ * 4 * votca::tools::conv::Pi) *
        dist_.data().y().cwiseQuotient(dist_.data().x().cwiseAbs2());

  } else {
    dist_.data().y() = scale_ /
                       ((double)nframes * area_ * rmax_ / (double)nbin_) *
                       dist_.data().y();
  }
  dist_.data().Save(out_ + suffix);
}

namespace Eigen {
std::istream &operator>>(std::istream &in, Vector3d &v) {
  char c;
  in.get(c);
  if (c != '[') {
    throw std::runtime_error("error, invalid character in vector string");
  }

  std::string str;
  while (in.good()) {
    in.get(c);
    if (c == ']') {  // found end of vector
      std::vector<double> d =
          votca::tools::Tokenizer(str, ",").ToVector<double>();
      if (d.size() != 3) {
        throw std::runtime_error("error, invalid number of entries in vector");
      }
      v.x() = d[0];
      v.y() = d[1];
      v.z() = d[2];
      return in;
    }
    str += c;
  }
  throw std::runtime_error(
      "did not find closing bracket in string to vec conversion");

  return in;
}
}  // namespace Eigen

void CsgDensityApp::EndEvaluate() {
  if (block_length_ == 0) {
    WriteDensity(frames_);
  }
}

// add our user program options
void CsgDensityApp::Initialize() {
  CsgApplication::Initialize();
  // add program option to pick molecule
  AddProgramOptions("Specific options:")(
      "type",
      boost::program_options::value<string>(&dens_type_)->default_value("mass"),
      "density type: mass or number")(
      "axis",
      boost::program_options::value<string>(&axisname_)->default_value("r"),
      "[x|y|z|r] density axis (r=spherical)")(
      "step",
      boost::program_options::value<double>(&step_)->default_value(0.01),
      "spacing of density")("block-length",
                            boost::program_options::value<votca::Index>(),
                            "  write blocks of this length, the averages are "
                            "cleared after every write")(
      "out", boost::program_options::value<string>(&out_), "Output file")(
      "rmax", boost::program_options::value<double>(),
      "rmax (default for [r] =min of all box vectors/2, else l )")(
      "scale",
      boost::program_options::value<double>(&scale_)->default_value(1.0),
      "scale factor for the density")(
      "molname",
      boost::program_options::value<string>(&molname_)->default_value("*"),
      "molname")(
      "filter",
      boost::program_options::value<string>(&filter_)->default_value("*"),
      "filter bead names")(
      "ref", boost::program_options::value<Eigen::Vector3d>(&ref_),
      "reference zero point");
}
