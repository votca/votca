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

#include <boost/format.hpp>
#include <fstream>
#include <iostream>
#include <votca/csg/csgapplication.h>

using namespace votca::csg;
using namespace std;
using boost::format;

class GmxTopolApp : public CsgApplication {
public:
  string ProgramName() { return "csg_gmxtopol"; }
  void HelpText(ostream &out) {
    out << "Create skeleton for gromacs topology based on atomistic topology\n"
           "and a mapping file. File still needs to be modified by the user.";
  }

  bool DoMapping(void) { return true; }

  void Initialize(void);
  bool EvaluateOptions(void) {
    CsgApplication::EvaluateOptions();
    CheckRequired("out", "no output topology specified");
    return true;
  }
  bool EvaluateTopology(Topology *top, Topology *top_ref);

protected:
  void WriteAtoms(ostream &out, Molecule &cg);
  void WriteInteractions(ostream &out, Molecule &cg);
  void WriteMolecule(ostream &out, Molecule &cg);
};

void GmxTopolApp::Initialize(void) {
  CsgApplication::Initialize();
  AddProgramOptions()(
      "out", boost::program_options::value<string>(),
      "output topology (will create .top and in future also .itp)");
}

bool GmxTopolApp::EvaluateTopology(Topology *top, Topology *top_ref) {
  if (top->MoleculeCount() > 1)
    cout << "WARNING: cannot create topology for topology with"
            "multiple molecules, using only first molecule\n";
  ofstream fl;
  fl.open((OptionsMap()["out"].as<string>() + ".top").c_str());
  WriteMolecule(fl, *(top->MoleculeByIndex(0)));
  fl.close();
  return true;
}

void GmxTopolApp::WriteAtoms(ostream &out, Molecule &cg) {
  out << "[atoms]\n";
  out << "; nr type resnr residue atom cgnr charge mass\n";
  for (int i = 0; i < cg.BeadCount(); ++i) {
    Bead *b = cg.getBead(i);
    out << format("%d %s 1 RES %s %d %f %f\n") % (i + 1) % b->getType() %
               b->getName() % (i + 1) % b->getQ() % b->getMass();
  }
  out << endl;
}

void GmxTopolApp::WriteInteractions(ostream &out, Molecule &cg) {
  int nb = -1;

  Interaction *ic;
  vector<Interaction *>::iterator iter;

  InteractionContainer &ics = cg.getParent()->BondedInteractions();

  for (iter = ics.begin(); iter != ics.end(); ++iter) {
    ic = *iter;
    if (ic->getMolecule() != cg.getId())
      continue;
    if (nb != ic->BeadCount()) {
      nb = ic->BeadCount();
      switch (nb) {
      case 2:
        out << "\n[ bonds ]\n";
        break;
      case 3:
        out << "\n[ angles ]\n";
        break;
      case 4:
        out << "\n[ dihedrals ]\n";
        break;
      default:
        string err = "cannot handle number of beads in interaction:";
        err += to_string(ic->getMolecule() + 1) + ":" + ic->getGroup();
        err += ":" + to_string(ic->getIndex() + 1);
        throw runtime_error(err);
      }
    }
    for (int i = 0; i < nb; ++i)
      out << ic->getBeadId(i) + 1 << " ";
    out << "  1  ; ";
    out << to_string(ic->getMolecule() + 1);
    out << ":" + ic->getGroup();
    out << ":" + to_string(ic->getIndex() + 1) << endl;
  }
}

void GmxTopolApp::WriteMolecule(ostream &out, Molecule &cg) {
  out << "[ moleculetype ]\n";
  out << cg.getName() << " 3\n\n";

  WriteAtoms(out, cg);
  WriteInteractions(out, cg);
}

int main(int argc, char **argv) {
  GmxTopolApp app;
  return app.Exec(argc, argv);
}
