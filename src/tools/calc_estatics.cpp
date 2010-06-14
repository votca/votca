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
        cout<<"TEstatic energy for neutral charged and diff: "<<crged<<" "<<neutr<<" "<<crged - neutr<<"\n";
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
        if(*imol == mol) continue;
        for(int i=0; i<mol->BeadCount();i++)
                for(int j=0; j!= (*imol)->BeadCount();j++){
            Bead *bi =  mol->getBead(i); Bead *bj = (*imol)->getBead(j);
            //vec r_v = atop->BCShortestConnection(bi->getPos(), bj->getPos());
            vec r_v=bi->getPos()-bj->getPos();
            double r=abs(r_v);
            double qi=bi->getQ(); double qj=bj->getQ();
            //cout<<bj->getPos()<<r_v<<r<<qi<<qj<<epsilon_dielectric<<pot<<qi*qj/epsilon_dielectric*1/r<<endl;
            pot+=qi*qj/epsilon_dielectric*1/r;
            //cout<<pot<<endl;
        }
    }
    return pot;
}


//Calculate Estatics Sasha
/*double CalcEstatics::CalcPot2(Topology *atop, Molecule *mol)
{
    double epsilon_dielectric=1.0;
    double pot=0.0;
    
    MoleculeContainer::iterator imol;
    for(imol = atop->Molecules().begin();imol!=atop->Molecules().end();imol++) {
        if(*imol == mol) continue;
	// charged
        CrgUnit *crg1 = mol->getUserData<CrgUnit>();
	// neutral
        CrgUnit *crg2 = (*imol)->getUserData<CrgUnit>();
	vec r_1 = crg1->GetCom();
	vec oldcom = crg2->GetCom();
	vec r_12 = atop->BCShortestConnection(r_1, oldcom); 
	double dist_12 = sqrt(r_12.x()*r_12.x()+r_12.y()*r_12.y()+r_12.z()*r_12.z());

	vector <vec> oldpos;
	vector <vec> ::iterator ipos;
	int count = 0;
	for(ipos=crg2->GetPosFirst(); ipos!=crg2->GetPosLast(); ++ipos, ++count){
	    oldpos.push_back(*ipos);
	    vec newpos = r_12 + crg2->GetPos(count) + crg1->GetCom() - oldcom;
	    crg2->SetPos(count, newpos);
	}

	for(int i=0; i < crg1->GetN(); i++) {
		for(int j=0; j < crg2->GetN(); j++) {
			vec r_v = atop->BCShortestConnection(crg1->GetPos(i), crg2->GetPos(j));
			double r=abs(r_v);
			double qi = mol->getBead(i)->getQ();	//  Error-prone! Change this!
			double qj = (*imol)-getBead(j)->getQ(); //  Error-prone! Change this!
			pot+=qi*qj/epsilon_dielectric*1/r;
		}
	}

	// Copy back
	count = 0;
	for(ipos=crg2->GetPosFirst(); ipos!=crg2->GetPosLast(); ++ipos, ++count){
	    crg2->SetPos(count, oldpos[count]);
	}
    }
    return pot;
}*/
