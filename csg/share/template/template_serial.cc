/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#include <cstdlib>
#include <votca/csg/beadlist.h>
#include <votca/csg/csgapplication.h>
#include <votca/csg/nblist.h>
#include <votca/csg/nblistgrid.h>
#include <votca/tools/histogram.h>

// using namespace votca::tools;
using namespace std;
using namespace votca::csg;

class CsgTestApp : public CsgApplication {
  string ProgramName() override { return "template_nblist"; }
  void HelpText(ostream &out) override {
    out << "rough template for rdf calculations";
  }

  void Initialize() override;

  bool DoTrajectory() override { return true; }

  void BeginEvaluate(Topology *top, Topology *top_ref) override;
  void EvalConfiguration(Topology *top, Topology *top_ref) override;
  void EndEvaluate() override;

 protected:
  votca::tools::HistogramNew rdf_;
  double cut_off_;
};

int main(int argc, char **argv) {
  CsgTestApp app;

  return app.Exec(argc, argv);
}

void CsgTestApp::EvalConfiguration(Topology *top, Topology *) {
  BeadList b;
  b.Generate(*top, "*");
  NBListGrid nb;
  nb.setCutoff(cut_off_);
  nb.Generate(b);
  for (auto &pair : nb) {
    rdf_.Process(pair->dist());
  }
}

void CsgTestApp::Initialize() {
  CsgApplication::Initialize();
  AddProgramOptions("RDF options")(
      "c", boost::program_options::value<double>()->default_value(1.0),
      "the cutoff");
}

void CsgTestApp::BeginEvaluate(Topology *, Topology *) {
  cut_off_ = OptionsMap()["c"].as<double>();
  rdf_.Initialize(0, cut_off_, 50);
}

void CsgTestApp::EndEvaluate() { rdf_.data().Save("rdf.dat"); }
