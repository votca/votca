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
        out << "Map a reference trajectory to a coarse-grained trajectory.\n"
            "This program can be used to map a whole trajectory or to\n"
            "create an initial configuration for a coarse-grained run only.";
    }

    bool DoTrajectory() { return true;}
    bool DoMapping() { return true;}

    void Initialize() {
        CsgApplication::Initialize();
        AddProgramOptions()
            ("out", boost::program_options::value<string>(),
                "  output file for coarse-grained trajectory")
                ("hybrid", "Create hybrid trajectory containing both atomistic and coarse-grained");
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
            // simply write the topology mapped by csgapplication classe
            _writer->Write(top);
        } else {
            // we want to combinge atomistic and coarse-grained into one topology
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
    string _hybrid_string;

};

void CsgMapApp::BeginEvaluate(Topology *top, Topology *top_atom) {
    string out = OptionsMap()["out"].as<string > ();
    string hybrid = OptionsMap()["hybrid"].as<string > ();
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


    _writer->Open(out);
};

int main(int argc, char **argv)
{
    CsgMapApp app;
    return app.Exec(argc, argv);
}

