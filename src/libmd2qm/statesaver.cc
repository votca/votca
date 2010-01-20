
#include <stdio.h>
#include <string>
#include "statesaver.h"

using namespace std;

StateSaver::StateSaver(QMTopology &qmtop) {
    _qmtop = &qmtop;
}

//void StateSaver::Open(string file)

void StateSaver::Save(string file, bool bAppend) {
    //_out = fopen(file.c_str(), bAppend ? "at" : "wt");
    _out.open(file.c_str(), ios::out | ios::binary);
    Write_QMBeads();
    //Write NBL
    _out.close();
    //_out = datafile(datafile.c_str(), ios::binary|ios::out);
}

void StateSaver::Write_QMBeads() {
    assert(_out.is_open());

    write<unsigned long>(_qmtop->BeadCount());
    for (BeadContainer::iterator iter = _qmtop->Beads().begin();
            iter != _qmtop->Beads().end(); ++iter) {
        QMBead *bi = dynamic_cast<QMBead*> (*iter);

        write<byte_t > (bi->getSymmetry());
        write<string > (bi->getName());
        write<string > (bi->getType()->getName());
        write<int>(bi->getResnr());
        write<double>(bi->getM());
        write<double>(bi->getQ());

        write<string > (bi->GetCrgUnit()->GetType()->GetName());
        write<unsigned short>(bi->getiPos());
        write<vec > (bi->getPos());
        write<vec > (bi->getU());
        write<vec > (bi->getV());
    }


}

void StateSaver::Load(string file) {
    _qmtop->Cleanup();
    _in.open(file.c_str(), ios::in | ios::binary);
    Read_QMBeads();
    //Read NBL
    _in.close();
}

void StateSaver::Read_QMBeads() {
    assert(_in.is_open());
   
    unsigned long nr_qmbeads = read<unsigned long>();
    cout << "Total number of QMBeads is " << nr_qmbeads << "\n";
    for (unsigned long i = 0; i < 1; i++) {

        byte_t symmetry = read<byte_t> ();
        string bead_name = read<string> ();
        string type_name = read<string> ();
        BeadType *type = _qmtop->GetOrCreateBeadType(type_name);
        int resnr = read<int>();
        double M = read<double>();
        double Q = read<double>();

        string crg_unit_name = read<string> ();
        unsigned short ipos = read<unsigned short>();
        vec Pos = read<vec> ();
        vec U = read<vec> ();
        vec V = read<vec> ();

        cout << "Bead Symmetry " << symmetry << "\n";
        cout << "Bead Name " << bead_name << "\n";
        cout << "Bead Type " << type_name << "\n";
        cout << "Residue Number " << resnr << "\n";
        cout << "Bead Mass " << M << "\n";
        cout << "Bead Charge " << Q << "\n";

        Bead *bead = _qmtop->CreateBead(symmetry, bead_name, type, resnr, M, Q);
        int molid = bead->getMolecule()->getId();
        string molandtype = lexical_cast<string > (molid) + ":" + crg_unit_name;
        cout << "molandtype " << molandtype << "\n";
        // CrgUnit * acrg = GetCrgUnitByName(molandtype);
        // if(acrg == NULL)
        //     acrg = CreateCrgUnit(molandtype, crg_unit_name, molid);

        // bead->setCrg(acrg);
        // bead->setPos(ipos);
        bead->setPos(Pos);
        bead->setU(U);
        bead->setV(V);

        cout << "The charge unit is called " << crg_unit_name << "\n";
        cout << "This bead is at int position " << ipos << "\n";
        cout << "This bead hast position  " << U << "\n";
        cout << "This bead has U " << U << "\n";
        cout << "This bead has V " << V << "\n";
                _in.close();
    }
}

