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

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <votca/csg/csgapplication.h>
#include <votca/csg/trajectorywriter.h>

using namespace std;
using namespace votca::csg;

class CsgMapApp
    : public CsgApplication
{
public:
    string ProgramName() { return "csg_map"; }
    void HelpText(ostream &out) {
      out << "Convert a reference atomistic trajectory or configuration into a coarse-grained one \n"
          << "based on a mapping xml-file. The mapping can be applied to either an entire trajectory \n"
          << "or a selected set of frames only (see options). The conventional file extensions of \n"
	  << "the supported simulation engines, e.g. Gromacs, are recognized. In the case of DL_POLY \n"
          << "the following convention is adopted: <name>.dlpf, <name>.dlpc and <name>.dlph  must \n"
	  << "be used for topology, configuration and trajectory files (aka FIELD, CONFIG and HISTORY). \n"
	  << "Moreover, by convention, <name> can be omitted, i.e. stand-alone extensions can be used, \n"
	  << "which invoke the corresponding standard names: FIELD, CONFIG, HISTORY (in the input), \n"
	  << "or FIELD_CGV, CONFIG_CGV, HISTORY_CGV (in the output). \n"
	  << "- NOTE: the tool can also be used for cross-converting configuration or trajectory files \n"
	  << "between the formats of the supported simulation engines!\n\n"
	  << "Examples:\n"
	  << "* csg_map --top FA-gromacs.tpr --trj FA-gromacs.trr --out CG-gromacs.xtc --cg cg-map.xml\n"
	  << "  - convert FA-gromacs.trr to CG-gromacs.xtc (trajectory in Gromacs format)\n"
	  << "* csg_map --top FA-gromacs.tpr --trj FA-gromacs.gro --out CG-gromacs.gro --cg cg-map.xml\n"
	  << "  - convert FA-gromacs.gro to CG-gromacs.gro (configuration in Gromacs format)\n"
	  << "* csg_map --top FA-dlpoly.dlpf --trj FA-dlpoly.dlph --out CG-dlpoly.dlph --cg cg-map.xml\n"
	  << "  - convert FA-dlpoly.dlph to CG-dlpoly.dlph (trajectory in DL_POLY format)\n"
	  << "* csg_map --top FA-gromacs.tpr --trj FA-gromacs.gro --out CG-dlpoly.dlpc --no-map\n"
	  << "  - convert FA-gromacs.gro to CG-dlpoly.dlpc (configuration: Gromacs to DL_POLY format)\n"
	  << "* csg_map --top .dlpf --trj .dlph --out .dlph --cg cg-map.xml --nframes 100\n"
	  << "  - convert HISTORY to HISTORY_CGV (trajectory in DL_POLY format: 100 frames only)\n";
    }

    bool DoTrajectory() { return true;}
    bool DoMapping() { return true;}

    void Initialize() {
        CsgApplication::Initialize();
        AddProgramOptions()
	  ("out", boost::program_options::value<string>(), 
           "  output file for coarse-grained trajectory; examples:\n"
           "  - <name>.trr, <name>.xtc, <name>.gro (Gromacs),\n"
           "  - <name>.dlph, <name>.dlpc, .dlph, .dlpc (DL_POLY)\n  ('.dlph'='use HISTORY_CGV', '.dlpc'='use CONFIG_CGV')")
	  ("vel", "  Write mapped velocities (if available)")
	  ("hybrid", "  Create hybrid trajectory containing both atomistic and coarse-grained");
    }

    bool EvaluateOptions() {
        CsgApplication::EvaluateOptions();
        CheckRequired("trj", "no trajectory file specified");
        CheckRequired("out", "need to specify output trajectory");
        return true;
    }

    void BeginEvaluate(Topology *top, Topology *top_ref);
void EvalConfiguration(Topology *top, Topology *top_ref) {
        if (!_do_hybrid) {
            // simply write the topology mapped by csgapplication class
            if (_do_vel) top->SetHasVel(true);
            _writer->Write(top);
        } else {
            // we want to combine atomistic and coarse-grained into one topology
            Topology *hybtol = new Topology();

            BeadContainer::iterator it_bead;
            ResidueContainer::iterator it_res;
            MoleculeContainer::iterator it_mol;
            InteractionContainer::iterator it_ia;

            hybtol->setBox(top->getBox());
            hybtol->setTime(top->getTime());
            hybtol->setStep(top->getStep());

            // copy all residues from both
            for (it_res = top_ref->Residues().begin(); it_res != top_ref->Residues().end(); ++it_res) {
                hybtol->CreateResidue((*it_res)->getName());
            }
            for (it_res = top->Residues().begin(); it_res != top->Residues().end(); ++it_res) {
                hybtol->CreateResidue((*it_res)->getName());
            }

            // copy all molecules and beads
          
            for(it_mol=top_ref->Molecules().begin();it_mol!=top_ref->Molecules().end(); ++it_mol) {
                Molecule *mi = hybtol->CreateMolecule((*it_mol)->getName());
                for (int i = 0; i < (*it_mol)->BeadCount(); i++) {
                    // copy atomistic beads of molecule
                    int beadid = (*it_mol)->getBead(i)->getId();

                    Bead *bi = (*it_mol)->getBead(i);
                    BeadType *type = hybtol->GetOrCreateBeadType(bi->getType()->getName());
                    Bead *bn = hybtol->CreateBead(bi->getSymmetry(), bi->getName(), type, bi->getResnr(), bi->getM(), bi->getQ());
                    bn->setOptions(bi->Options());
                    bn->setPos(bi->getPos());
                    if (bi->HasVel()) bn->setVel(bi->getVel());

                    mi->AddBead(hybtol->Beads()[beadid], (*it_mol)->getBeadName(i));

                }

                if (mi->getId() < top->MoleculeCount()) {
                    // copy cg beads of molecule
                    Molecule *cgmol = top->Molecules()[mi->getId()];
                    for (int i = 0; i < cgmol->BeadCount(); i++) {
                        Bead *bi = cgmol->getBead(i);
                        // todo: this is a bit dirty as a cg bead will always have the resid of its first parent
                        Bead *bparent = (*it_mol)->getBead(0);
                        BeadType *type = hybtol->GetOrCreateBeadType(bi->getType()->getName());
                        Bead *bn = hybtol->CreateBead(bi->getSymmetry(), bi->getName(), type, bparent->getResnr(), bi->getM(), bi->getQ());
                        bn->setOptions(bi->Options());
                        bn->setPos(bi->getPos());
                        if (bi->HasVel()) bn->setVel(bi->getVel());
                        int mid = bparent->getMolecule()->getId();
                        mi->AddBead(bi, bi->getName());
                    }
                }
                
            }
           hybtol->setBox(top_ref->getBox());

            _writer->Write(hybtol);
        }
    }

    void EndEvaluate() {
        _writer->Close();
        delete _writer;
    }

protected:
    TrajectoryWriter *_writer;
    bool _do_hybrid;
    bool _do_vel;

};

void CsgMapApp::BeginEvaluate(Topology *top, Topology *top_atom) {
    string out = OptionsMap()["out"].as<string > ();
    cout << "writing coarse-grained trajectory to " << out << endl;
    _writer = TrjWriterFactory().Create(out);
    if (_writer == NULL)
        throw runtime_error("output format not supported: " + out);

    _do_hybrid = false;
    if(OptionsMap().count("hybrid")){
        if (!_do_mapping)
            throw runtime_error("options hybrid and no-map not compatible");
        cout << "Doing hybrid mapping..." << endl;
        _do_hybrid = true;
    }

    _do_vel = false;
    if(OptionsMap().count("vel")){
        _do_vel = true;
    }

    _writer->Open(out);
};

int main(int argc, char **argv)
{
    CsgMapApp app;
    return app.Exec(argc, argv);
}

