#include "sandbox.h"

namespace votca { namespace ctp {

    void Sandbox::Initialize(QMTopology *top, Property *opt) {

        cout << "Initialize (Sandbox)..." << endl;

        string key;

        key = "options.sandbox.id";
        if (opt->exists(key)) {
            _ID = opt->get(key).as< int >();
        }

        key = "options.sandbox.sec1";
        if (opt->exists(key+".p1")) {
            _p1 = opt->get(key+".p1").as< double >();
        }

        key = "options.sandbox.sec2";
        if (opt->exists(key+".p2")) {
            _p2 = opt->get(key+".p2").as< double >();
        }

        cout << "P1     " << _p1 << endl;
        cout << "P2     " << _p2 << endl;
        cout << "ID     " << _ID << endl;

    }

    bool Sandbox::EvaluateFrame(QMTopology *top) {

        cout << "Calculate (Sandbox)... " << endl;

        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        cout << "Bead iterator" << endl;
        vector < Bead* >::iterator ait;
        for (ait = top->Beads().begin();
             ait != top->Beads().end();
             ait++) {
            Bead* atm = *ait;
            cout << "getID    " << atm->getId() << endl;
            cout << "getResnr " << atm->getResnr() << endl;
            cout << "getQ     " << atm->getQ() << endl;
        }

        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        cout << "Molecule iterator" << endl;
        vector < Molecule* >::iterator chuit;
        for (chuit = top->Molecules().begin();
             chuit != top->Molecules().end();
             chuit ++) {
            Molecule* mol = *chuit;
            cout << "getName    " << mol->getName() << endl;
            cout << "getID      " << mol->getId() << endl;
            cout << "Bead count " << mol->BeadCount() << endl;
        }

        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        cout << "Charge-unit iterator" << endl;
        vector < QMCrgUnit* >::iterator chrgit;
        for (chrgit = top->CrgUnits().begin();
             chrgit != top->CrgUnits().end();
             chrgit++ ) {
            QMCrgUnit* unit = *chrgit;
            cout << "getType->GetName " << (unit->getType())->GetName() << endl;
            cout << "getID            " << unit->getId() << endl;
            cout << "getName          " << unit->getName() << endl;
        }

        vector < QMCrgUnit* > lcharges = top->CrgUnits();
        Topology atop;
        atop.setBox(top->getBox());
        for (vector<QMCrgUnit*>::iterator itl = lcharges.begin();
                itl != lcharges.end();
                itl++ ) {
            top->AddAtomisticBeads(*itl, &atop);
        }
        cout << "# chrg units in QMTOP " << lcharges.size() << endl;
        cout << "# molecules in ATMTOP " << atop.MoleculeCount() << endl;

        MoleculeContainer::iterator molit;
        for (molit = atop.Molecules().begin();
                molit != atop.Molecules().end();
                molit++ ) {
            Molecule *mol = *molit;
            cout << "# BeadCount() in " << mol->getName();
            cout << ": " << mol->BeadCount() << endl;
        
            for (int i = 0;
                     i < mol->BeadCount();
                     i++ ) {
                
                Bead *bead = mol->getBead(i);
                cout << " Bead ID   " << bead->getId();
                cout << " Bead Pos  " << bead->getPos();
                cout << " Bead Q    " << bead->getQ();
                cout << " Bead Name " << bead->getName() << endl;
            }
        }


    }



}} /* exit namespace votca::ctp */
