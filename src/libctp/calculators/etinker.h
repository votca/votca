#ifndef _ETINKER_H
#define	_ETINKER_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>

#include <stdlib.h>
#include <math.h>
#include <list>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>

namespace votca { namespace ctp {

/**
    \brief Tinker input: xyz coordinates [Angstroem] with a given molecule centered in the box.

Callname: etinker

Part of the input of the TINKER program which is used to evaluate the polarization contribution to site energies (self-consistently). Dumps the coordinates [xyz, Angstroem] of all atoms in the snapshot or atom within a cutoff (nm, default 50nm) based on centers of mass of molecules. Files are named xyz_N, where N (starts at 0, first molecule in the file) is the number of the molecule whose site energy is computed. This molecule is placed in the center of the box and the nearest image convention is used for the rest of molecules.

*/

class Etinker : public QMCalculator
{
public:
    Etinker() {};
    ~Etinker() {};

    const char *Description() { return "Tinker input: xyz coordinates [Angstroem] with a given molecule centered in the box."; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void WriteAtoms(Topology *atop, Molecule *mol);
   
private:
    double _dump_cutoff;
    Property * _options;
};

inline void Etinker::Initialize(QMTopology *top, Property *options) {
    _options = options;
    if (options->exists("options.tinker.cutoff")) {
        _dump_cutoff = options->get("options.tinker.cutoff").as<double>();
        cout << "Writing out atomic XYZ coordinates for molecular cutoff: " << _dump_cutoff <<" nm"<<endl;
    } else {
        _dump_cutoff=50.0;
        cout << "Warning: No cutoff for molecules has been provided, using default cutoff of 50nm." << endl;
}

}

inline bool Etinker::EvaluateFrame(QMTopology *top) {
    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    Topology atop;
    atop.setBox(top->getBox());
    for (vector<QMCrgUnit *>::iterator itl = lcharges.begin(); itl != lcharges.end(); itl++) {
        top->AddAtomisticBeads(*itl, &atop);
    }
    cout << "Number of charge units in top " << lcharges.size() << endl;
    cout << "Number of molecules in atop " << atop.MoleculeCount() << endl;
    double crged, neutr;

    MoleculeContainer::iterator imol;
    for (imol = atop.Molecules().begin(); imol != atop.Molecules().end(); imol++) {
        Molecule *mol = *imol;
        CrgUnit *crg = mol->getUserData<CrgUnit > ();
        WriteAtoms(&atop, mol);
        //cout << "Estatic energy [eV] for charged / neutral / crg-neutr=espilon: " << crged << " " << neutr << " " << crged - neutr << "\n";
        cout << "writing out coordinats in Angstroem for crgunit " << crg->getId()<<"\n";
    }
    return true;
}


inline void Etinker::WriteAtoms(Topology *atop, Molecule *mol) //wegen Übergabe per * unten ->
{
    MoleculeContainer::iterator cmol;
    int nr = mol->getId();
            string filename = "xyz_" + boost::lexical_cast<string > (nr);
            FILE * data;
            //First write out the xyz coordinates of atoms belonging to molecule N which is in the center of the box in file xyz_N
            data=fopen(filename.c_str(),"w");
            for (int j = 0; j != mol->BeadCount(); j++) {
                Bead *bj = mol->getBead(j);
                vec r_v = mol->getUserData<CrgUnit > ()->GetCom()-(bj->getPos());
                fprintf(data, "%d %d %.6f %.6f %.6f\n",mol->getId(),mol->getId(),r_v.getX()*10.0,r_v.getY()*10.0,r_v.getZ()*10.0);
        }
     MoleculeContainer::iterator dmol;
         for (dmol = atop->Molecules().begin(); dmol != atop->Molecules().end(); dmol++) {
           if (*dmol == mol) continue;
        //We should replace getUserData by a funtion GetCom for molecules
        vec bcs = atop->BCShortestConnection(mol->getUserData<CrgUnit > ()->GetCom(), (*dmol)->getUserData<CrgUnit > ()->GetCom());
         if ( abs(bcs)> _dump_cutoff) continue;

        vec dist = (*dmol)->getUserData<CrgUnit > ()->GetCom() - mol->getUserData<CrgUnit > ()->GetCom();
        vec diff = bcs - dist;
        for (int j = 0; j != (*dmol)->BeadCount(); j++) {
                Bead *bj = (*dmol)->getBead(j);
                //distance vector may have to be shifted by s
                vec r_v = mol->getUserData<CrgUnit > ()->GetCom()-(bj->getPos() + diff);
               // data << mol->getId()<<" "<<(*dmol)->getId()<<" "<< r_v.getX()*10.0 <<" "<< r_v.getY()*10.0 <<" "<< r_v.getZ()*10.0 <<"\n";
                 fprintf(data, "%d %d %.6f %.6f %.6f\n",mol->getId(),(*dmol)->getId(),r_v.getX()*10.0,r_v.getY()*10.0,r_v.getZ()*10.0);
        }
        
    }
    fclose(data);
    //data.close();
}

}}

#endif	/* _ETINKER_H */

