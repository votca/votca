#include <stdlib.h>
#include "dump_atoms.h"
#include <math.h>
#include <list>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>

void DumpAtoms::Initialize(QMTopology *top, Property *options) {
    _options = options;
    if (options->exists("options.dump_atoms_params.cutoff")) {
        _dump_cutoff = options->get("options.dump_atoms_params.cutoff").as<double>();
        cout << "Writing out atomic XYZ coordinates for molecular cutoff: " << _dump_cutoff <<" nm"<<endl;
    } else {
        _dump_cutoff=50.0;
        cout << "Warning: No cutoff for molecules has been provided, using default cutoff of 50nm." << endl;
}

}

bool DumpAtoms::EvaluateFrame(QMTopology *top) {
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


void DumpAtoms::WriteAtoms(Topology *atop, Molecule *mol) //wegen Übergabe per * unten ->
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

