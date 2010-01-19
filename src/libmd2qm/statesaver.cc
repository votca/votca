
#include <stdio.h>
#include <string>
#include "statesaver.h"

using namespace std;
StateSaver::StateSaver(QMTopology& qmtop)
{}

//void StateSaver::Open(string file)
void StateSaver::Open(string file, bool bAppend)
{
    //_out = fopen(file.c_str(), bAppend ? "at" : "wt");
    ofstream _out(file.c_str(), ios::binary|ios::out);
    //_out = datafile(datafile.c_str(), ios::binary|ios::out);
}

void StateSaver::Close()
{
 _out.close();
}

void StateSaver::Write_QMBeads(QMTopology *top)
{   QMTopology *qmtop = top;
    for(BeadContainer::iterator iter=qmtop->Beads().begin();
    iter!=qmtop->Beads().end(); ++iter) {
        QMBead *bi = dynamic_cast<QMBead*>(*iter);

    string crg_unit_name = bi->GetCrgUnit()->GetType()->GetName();
        unsigned short len_name = crg_unit_name.length();
    int ipos=bi->getiPos();
    string type=bi->getType()->getName();
        unsigned short len_type = type.length();
    int symmetry=bi->getSymmetry();
    string bead_name=bi->getName();
        unsigned short len_bead_name = bead_name.length();
    int res_number=bi->getResnr();
    double mass=bi->getM();
    double charge=bi->getQ();
    vec r = bi->getPos();
    vec u = bi->getU();
    vec v = bi->getV();

    //writing all out
     _out.write((char *) &len_name, sizeof(len_name));
     cout<<"We have"<<len_name<<"letters in the name\n";
     _out.write((char *) crg_unit_name.c_str(), sizeof(crg_unit_name));
     _out.write((char *) &ipos, sizeof(ipos));
     _out.write((char *) &len_type, sizeof(len_type));
     _out.write((char *) type.c_str(), sizeof(type));
     _out.write((char *) &symmetry, sizeof(symmetry));
     _out.write((char *) &len_bead_name, sizeof(len_bead_name));
     _out.write((char *) bead_name.c_str(), sizeof(bead_name));
     _out.write((char *) &res_number, sizeof(len_bead_name));
     _out.write((char *) &mass, sizeof(mass));
     _out.write((char *) &charge, sizeof(charge));
     _out.write((char *) &r, sizeof(r));
     _out.write((char *) &u, sizeof(u));
     _out.write((char *) &v, sizeof(v));
    }

    // fflush(_out);
}

/*VICTOR:
void PDBWriter::Open(string file, bool bAppend)
{
    _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void PDBWriter::Close()
{
    fclose(_out);
}

void PDBWriter::Write(Topology *conf)
{
    Topology *top = conf;
    fprintf(_out, "MODEL     %4d\n", conf->getStep());
    for(BeadContainer::iterator iter=conf->Beads().begin();
    iter!=conf->Beads().end(); ++iter) {
        Bead *bi = *iter;
        vec r = bi->getPos();
        fprintf(_out,
                "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                (bi->getId()+1)%100000,   // atom serial number
                bi->getName().c_str(),  // atom name
                top->getResidue(bi->getResnr())->getName().c_str(), // residue name
                " ", // chain identifier 1 char
                bi->getResnr()+1, // residue sequence number
                10.*r.x(), 10.*r.y(), 10.*r.z(),
                bi->getQ(), bi->getM());  // is this correct??

        if(bi->getSymmetry()>=2) {
           vec ru = 0.1*bi->getU() + r;

            fprintf(_out,
                "HETATM%5d %4s %3s %1s%4d    %8.3f%8.3f%8.4f%6.2f%6.2f\n",
                bi->getId()+1,   // atom serial number
                bi->getName().c_str(),  // atom name
                "REU", // residue name
                " ", // chain identifier 1 char
                bi->getResnr()+1, // residue sequence number
                10.*ru.x(), 10.*ru.y(), 10.*ru.z(),
                0., 0.);  // is this correct??
        }
        if(bi->getSymmetry()>=3) {
           vec rv = 0.1*bi->getV() + r;
            fprintf(_out,
                "HETATM%5d %4s %3s %1s%4d    %8.3f%8.3f%8.4f%6.2f%6.2f\n",
                bi->getId()+1,   // atom serial number
                bi->getName().c_str(),  // atom name
                "REV", // residue name
                " ", // chain identifier 1 char
                bi->getResnr()+1, // residue sequence number
                10.*rv.x(), 10.*rv.y(), 10.*rv.z(),
                0.,0.);
        }
   }
    fprintf(_out, "ENDMDL\n");
    fflush(_out);
}*/
