/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

/**
    \brief class for writing dlpoly topology files

    This class encapsulates the dlpoly topology writing functions

*/

class DLPTopolApp
    : public CsgApplication
{
public:
    string ProgramName() { return "csg_dlptopol"; }
    void HelpText(ostream &out) {
        out << "Create a dlpoly topology template based on an existing (atomistic) topology and \n"
	    << "a mapping xml-file. The created template file needs to be inspected and amended by the user!\n\n"
	    << "Examples:\n"
	    << "* csg_dlptopol --top .dlpf --out .dlpf --cg cg-map.xml\n  convert FIELD to FIELD_CGV using cg-map.xml\n"
	    << "* csg_dlptopol --top FA-dlpoly.dlpf --out CG-dlpoly.dlpf --cg cg-map.xml\n"
	    << "* csg_dlptopol --top FA-gromacs.tpr --out FA-dlpoly.dlpf --no-map\n";
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
    void WriteMoleculeAtoms(ostream &out, Molecule &cg);
    void WriteMoleculeInteractions(ostream &out, Molecule &cg);
    void WriteVDWInteractions(ostream &out, Molecule &cg);
    void WriteMolecularType(ostream &out, Molecule &cg, int nummol);
};

void DLPTopolApp::Initialize(void)
{
    CsgApplication::Initialize();
    //AddProgramOptions()
    //("top", boost::program_options::value<string>(), 
    //"  input topology in any known format:\n  <name>.dlpf for dlpoly, <name>.tpr for gromacs\n  (convention: '.dlpf'='use FIELD')");
    AddProgramOptions()
      ("out", boost::program_options::value<string>(), 
       "  output topology in dlpoly format");
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

  // do CG mapping

  MoleculeContainer &mols = top->Molecules();
  MoleculeContainer MolecularTypes;
  MoleculeContainer::iterator iter;

  int prv_mol_number = 1;
  string prv_mol_name;
  vector<int> nummols;

  vector<string> vdw_pairs;

  // find all unique molecular types

  for(iter=mols.begin(); iter!=mols.end(); ++iter) {
    Molecule *mol = *iter;

    // molecules are ignored during the mapping stage 
    // i.e. the ignored ones do not enter the CG topology (*top) - ?
    //if( IsIgnored(mol->getName()) ) continue;

    if( mol->getName()==prv_mol_name ) {
      prv_mol_number++;
      continue;
    }

    nummols.push_back(prv_mol_number);
    prv_mol_number = 1;
    prv_mol_name   = mol->getName();

    //#ifdef DEBUG
    cout << "'" << mol->getName() << "' added to CG molecular types" << endl;
    //#endif

    MolecularTypes.push_back(mol);

    // collect unique bead pairs over all molecular types found

    for(int ib1=0; ib1<mol->BeadCount(); ib1++) {
      string bead_name1 mol->getBeadTypeName();
      bead_name1 = bead_name1.substr(0,bead_name1.find_first_of("#")); // skip #index of atom from its name

      for(unsigned int imt=0; imt<MolecularTypes.size(); imt++) {

	for(int ib2=0; ib2<MolecularTypes[imt]->BeadCount(); ib2++) {

    string bead_name2 = MolecularTypes[imt]->getBead(ib2)->getBeadTypeName();
    bead_name2 = bead_name2.substr(0,bead_name2.find_first_of("#")); // skip #index of atom from its name

	  stringstream ss_bp1,ss_bp2;

	  ss_bp1 << format("%8s%8s" ) % bead_name1 % bead_name2;
	  ss_bp2 << format("%8s%8s" ) % bead_name2 % bead_name1;

          bool is_new_pair=true;

	  for(unsigned int ibp=0; ibp<vdw_pairs.size(); ibp++) {
	    if( ss_bp1.str()==vdw_pairs[ibp] || ss_bp2.str()==vdw_pairs[ibp] ) { 
	      is_new_pair=false; 
	      break;
	    }
	  }
	  if( is_new_pair ) {
	    vdw_pairs.push_back(ss_bp1.str());
#ifdef DEBUG
	    cout << "'" << ss_bp1.str() << "' added to CG vdw pairs" << endl;
#endif
	  }
	}
      }
    }
  }
  nummols.push_back(prv_mol_number);

  if(MolecularTypes.size() > 1)
    cout << "WARNING: creation of topology for multiple molecular types "
      "is experimental at this stage\n";

  ofstream fl;
  fl.open(fname.c_str());

  fl << "From VOTCA with love" << " # check/amend this file if needed!\n";
  fl << "units kJ\n";
  fl << "molecular types " << MolecularTypes.size() << endl;

  for(unsigned int i=0; i<MolecularTypes.size(); i++) {
    WriteMolecularType(fl, *(MolecularTypes[i]), nummols[i+1]);
  }

  // vdw seciton (pairwise vdw/CG interactions)

  if( vdw_pairs.size()>0 ) {

    fl << "vdw "<< vdw_pairs.size() << endl;
   
    for(unsigned int ibp=0; ibp<vdw_pairs.size(); ibp++) {
      fl << vdw_pairs[ibp] << " tab   1.00000  0.00000\n";
    }
  }

  fl << "close" << endl;

  cout << "Created template for dlpoly topology - please, check & amend if needed!" << endl;

  fl.close();
  return true;
}


void DLPTopolApp::WriteMoleculeAtoms(ostream &out, Molecule &cg)
{
  out << "atoms " << cg.BeadCount() << endl;
    out << "# name  mass  charge  nrept  ifrozen (optional: ngroup, index, name/type, type/residue, index/res-ID) \n";
    for(int i=0; i<cg.BeadCount(); ++i) {
        Bead *b=cg.getBead(i);

        string bname=b->getName();
	      string btype = b->getBeadTypeName();

        bname = bname.substr(0,bname.find_first_of("#")); // skip #index of atom from its name
        btype = btype.substr(0,btype.find_first_of("#")); // skip #index of atom from its type

        out << format("%8s  %10f  %10f     1     0     1 %10d  %8s  %8s %10d \n")
            % btype % b->getMass() % b->getQ() % (i+1) % btype % bname % (i+1);
	//% b->getType()->getName() % b->getMass() % b->getQ() % (i+1) % b->getType()->getName() % b->getName() % (i+1);
    }
}

void DLPTopolApp::WriteMoleculeInteractions(ostream &out, Molecule &cg)
{
    InteractionContainer ics=cg.Interactions();
    vector<Interaction *>::iterator iter;

    stringstream sout;

    int n_entries = 0;
    int nb=-1;

    for(iter=ics.begin(); iter!=ics.end(); ++iter) {
        Interaction *ic = *iter;
        if(nb != ic->BeadCount()) {

	    if(sout.str()!="") 
	      out << n_entries << endl << sout.str();

	    sout.str("");
	    n_entries = 0;

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
                    string err = "cannot handle number of beads in interaction:";
                    err += to_string(ic->getMolecule()+1)+":"+ic->getGroup();
                    err += ":"+to_string(ic->getIndex()+1);
                    throw runtime_error(err);
            }
        }
	n_entries++;
	// to do: is it possible to use bond/angle/dihedral function types for 1:1 mapping? (CG map overwrites ic->Group anyway)
        //sout << ic->getInteractionFunc(); // something like that (only for 1:1 mapping!)
        sout << " tab ";
        for(int i=0; i<nb; ++i)
	  sout << ic->getBeadId(i)+1 << " ";
        sout << "   1.00000  0.00000" << " # ";
        sout << to_string(ic->getMolecule()+1) ;
        sout << ":"+ic->getGroup();
        sout << ":"+to_string(ic->getIndex()+1) << endl;
    }
    if(sout.str()!="") out << n_entries << endl << sout.str();
}

void DLPTopolApp::WriteMolecularType(ostream &out, Molecule &cg, int nummol)
{
    out << cg.getName() << endl;
    out << "nummols " << nummol << endl;

    WriteMoleculeAtoms(out, cg);
    WriteMoleculeInteractions(out, cg);

    out << "finish" << endl;
}

int main(int argc, char** argv)
{    
    DLPTopolApp app;
    return app.Exec(argc, argv);        
}

