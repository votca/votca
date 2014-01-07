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

#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include <votca/csg/csgapplication.h>

using namespace votca::csg;
using namespace std;
using boost::format;

class DLPTopolApp
    : public CsgApplication
{
public:
    string ProgramName() { return "csg_dlptopol"; }
    void HelpText(ostream &out) {
        out << "Create template for dlpoly topology (FIELD) based on atomistic one"
	  "and a mapping file (cg-map.xml).\n"
	  "File still needs to be inspected/modified by the user.";
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

void DLPTopolApp::Initialize(void)
{
    CsgApplication::Initialize();
    AddProgramOptions()
      ("out", boost::program_options::value<string>(), "output dlpoly topology (will create FIELD_CGV or <name>.dlpf )");
}

bool DLPTopolApp::EvaluateTopology(Topology *top, Topology *top_ref)
{
  // check the file names from options

  string fname=OptionsMap()["top"].as<string>();

  if( fname == ".dlpf" ) {
    fname = "FIELD";
  }

  cout << "input  file-name: " << fname << endl;

  fname=OptionsMap()["out"].as<string>();

#ifdef DEBUG
  cout << "output file-name given: " << fname << endl;
#endif

  if( fname == ".dlpf" ) {
    fname = "FIELD_CGV";
  }

#ifdef DEBUG
  cout << "output file-name actual: " << fname << endl;
#else
  cout << "output file-name: " << fname << endl;
#endif

  // do the mapping

  if(top->MoleculeCount() > 1)
    cout << "WARNING: cannot create topology for topology with"
      "multiple molecules, using only first molecule\n";
  ofstream fl;
  fl.open(fname.c_str());
  WriteMolecule(fl, *(top->MoleculeByIndex(0)));
  fl.close();
  return true;
}


void DLPTopolApp::WriteAtoms(ostream &out, Molecule &cg)
{
  out << "ATOMS " << cg.BeadCount() << "\n";
    out << "# name mass charge nrept ifrozen (optional: ngroup, index, type, residue) \n";
    for(int i=0; i<cg.BeadCount(); ++i) {
        Bead *b=cg.getBead(i);
       
        out << format("%4s %f %f 1 0 1 %d %s \n")
            % b->getType()->getName() % b->getM() % b->getQ() % (i+1) % b->getName();
    }
}

void DLPTopolApp::WriteInteractions(ostream &out, Molecule &cg)
{
    int nb=-1;
    
    Interaction *ic;
    vector<Interaction *>::iterator iter;
  
    InteractionContainer &ics=cg.getParent()->BondedInteractions();

    int n_entries = ics.end()-ics.begin()+1;

    for(iter=ics.begin(); iter!=ics.end(); ++iter) {
        ic = *iter;
        if(ic->getMolecule() != cg.getId()) continue;
        if(nb != ic->BeadCount()) {
            nb=ic->BeadCount();
            switch(nb) {
                case 2:
                    out << "bonds ";
                    break;
                case 3:
                    out << "angles ";
                    break;
                case 4:
                    out << "dihedrals ";
                    break;
                default:
                    throw runtime_error(string("cannot handle number of beads in interaction:") +
                            ic->getName());
            }
	    out << " -1  # check/amend the number of interactions!\n";
	    //out << n_entries << endl;
        }
        out << " tab ";
        for(int i=0; i<nb; ++i)
            out << ic->getBeadId(i)+1 << " ";
        out << " # " << ic->getName() << endl;
    }
}

void DLPTopolApp::WriteMolecule(ostream &out, Molecule &cg)
{
    out << "From VOTCA with love\n";
    out << "UNITS kJ\n";
    out << "MOLECULE TYPES 1\n";

    out << cg.getName() << endl;
    out << "NUMMOLS " << 1 << " # check/amend the number of molecule repetitions!\n";

    WriteAtoms(out, cg);
    WriteInteractions(out, cg);

    cout << "Created template for dlpoly topology - please, check & amend if needed!" << endl;
}

int main(int argc, char** argv)
{    
    DLPTopolApp app;
    return app.Exec(argc, argv);        
}

