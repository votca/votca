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
#include <csgapplication.h>
#include <trajectorywriter.h>

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
                ("hybrid", boost::program_options::value<string>(&_hybrid_string)->default_value(""),
                "  Create hybrid topology containing both atomistic and coarse-grained");
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

            // copy all residues
            for (it_res = top->Residues().begin(); it_res != top->Residues().end(); ++it_res) {
                hybtol->CreateResidue((*it_res)->getName());
                //cout << (*it_res)->getName() << endl;
            }
            // copy all molecules
            /// TODO: copy molecules below....
            for(it_mol=top->_molecules.begin();it_mol!=top->_molecules.end(); ++it_mol) {
                Molecule *mi = CreateMolecule((*it_mol)->getName());
                for(int i=0; i<(*it_mol)->BeadCount(); i++) {
                    int beadid = (*it_mol)->getBead(i)->getId();
                    mi->AddBead(_beads[beadid], (*it_mol)->getBeadName(i));
                }
            }

            // TODO: use bead-> getParent!

            BeadContainer::iterator iterAT = top_ref->Beads().begin();
            for (BeadContainer::iterator iterCG = top->Beads().begin();
                    (iterCG != top->Beads().end() && iterAT != top->Beads().end());
                    ++iterCG, ++iterAT) {
                Bead *beadAT = *iterAT;
                Bead *beadCG = *iterCG;
                //cout << beadAT->getName() << endl;
                //hybtol->CreateBead(1, "", top->GetOrCreateBeadType("no"), 0, 0, 0);
                cout << "AT " <<beadAT->getName() << beadAT->getPos() << endl;
                 Bead *hybBead= hybtol->CreateBead(beadAT->getSymmetry(), beadAT->getName(),
                 hybtol->GetOrCreateBeadType(beadAT->getType()->getName()), beadAT->getResnr(), beadAT->getM(),
                 beadAT->getQ());
                 hybBead->setPos(beadAT->getPos());
                 hybBead->setOptions(beadAT->Options());
                /*hybtol->CreateBead(beadCG->getSymmetry(), beadCG->getName(),
                beadCGet->Type(), beadAT->getResnr(), beadCG->getM(),
                beadCG->getQ());*/
            }
            for (it_bead = hybtol->Beads().begin(); it_bead != hybtol->Beads().end(); ++it_bead) {
                Bead *bi = *it_bead;
                cout << bi->getName() << bi->getPos() << endl;
            }
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

    if (hybrid == "yes") {
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

