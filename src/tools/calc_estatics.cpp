/*
 * File:   calc_estatics.cpp
 * Author: mayfalk
 *
 * Created on May 21, 2010, 11:21 PM
 */


#include <stdlib.h>
#include "calc_estatics.h"
#include <math.h>
#include <list>

void CalcEstatics::EvaluateFrame(QMTopology *top)
{   cout<<"Doing Estatics\n";
    list<CrgUnit *> lcharges = top->crglist();
    Topology atop;
    atop.setBox(top->getBox());
    for (list<CrgUnit *>::iterator itl = lcharges.begin(); itl != lcharges.end(); itl++){
        top->AddAtomisticBeads(*itl,&atop);
    }
    cout<<"Number of charge units in top "<<lcharges.size()<<endl;
    cout<<"Number of molecules in atop "<<atop.MoleculeCount()<<endl;
    MoleculeContainer::iterator imol;
    for(imol = atop.Molecules().begin();imol!=atop.Molecules().end();imol++) {
        Molecule *mol = *imol;
        
        //Compute the EstaticPotentialEnergy for each molecule mol if it is neutral
        double neutr = CalcPot(&atop, mol);

        CrgUnit *crg = mol->getUserData<CrgUnit>();

        top->CopyChargesOccupied(crg, mol);

        double crged = CalcPot(&atop,mol);
        crg->setEnergy(crged - neutr);

        top->CopyCharges(crg, mol);
        cout<<"Estatic energy for neutral charged and diff: "<<crged<<" "<<neutr<<" "<<crged - neutr<<"\n";
    }

    atop.Cleanup();

    
}

//Calculate Estatics
double CalcEstatics::CalcPot(Topology *atop, Molecule *mol) //wegen Übergabe per * unten ->
{
    double epsilon_dielectric=1.0;
    double pot=0.0;
    
    MoleculeContainer::iterator imol;
    for(imol = atop->Molecules().begin();imol!=atop->Molecules().end();imol++) {
        if(*imol != mol) continue;
        for(int i=0; i<mol->BeadCount();i++)
                for(int j=0; j!= (*imol)->BeadCount();j++){
            Bead *bi =  mol->getBead(i); Bead *bj = (*imol)->getBead(j);
            vec r_v = atop->BCShortestConnection(bi->getPos(), bj->getPos());
            //matrix box=atop->getBox();
            double r=abs(r_v);
            double qi=bi->getQ(); double qj=bj->getQ();
            cout<<bj->getPos()<<r_v<<r<<qi<<qj<<epsilon_dielectric<<pot<<qi*qj/epsilon_dielectric*1/r<<endl;
            pot+=qi*qj/epsilon_dielectric*1/r;
            cout<<pot<<endl;
        }
    }
    return pot;
    }


