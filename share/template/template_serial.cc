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

#include <stdlib.h>
#include <votca/csg/beadlist.h>
#include <votca/csg/csgapplication.h>
#include <votca/csg/nblist.h>
#include <votca/csg/nblistgrid.h>
#include <votca/tools/histogramnew.h>

// using namespace votca::tools;
using namespace std;
using namespace votca::csg;

class CsgTestApp : public CsgApplication {
  string ProgramName() { return "template_nblist"; }
  void HelpText(ostream &out) { out << "rough template for rdf calculations"; }

  void Initialize();

  bool DoTrajectory() { return true; }

  void BeginEvaluate(Topology *top, Topology *top_ref);
  void EvalConfiguration(Topology *top, Topology *top_ref);
  void EndEvaluate();

 protected:
  HistogramNew _rdf;
  double _cut_off;
};

int main(int argc, char **argv) {
  CsgTestApp app;

  return app.Exec(argc, argv);
}

void CsgTestApp::EvalConfiguration(Topology *top, Topology *top_ref) {
  BeadList b;
  b.Generate(*top, "*");
  NBListGrid nb;
  nb.setCutoff(_cut_off);
  nb.Generate(b);
  NBList::iterator i;
  for (i = nb.begin(); i != nb.end(); ++i) _rdf.Process((*i)->dist());
}

void CsgTestApp::Initialize() {
  CsgApplication::Initialize();
  AddProgramOptions("RDF options")(
      "c", boost::program_options::value<double>()->default_value(1.0),
      "the cutoff");
}

void CsgTestApp::BeginEvaluate(Topology *top, Topology *top_ref) {
  _cut_off = OptionsMap()["c"].as<double>();
  _rdf.Initialize(0, _cut_off, 50);
}

void CsgTestApp::EndEvaluate() { _rdf.data().Save("rdf.dat"); }
