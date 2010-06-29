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

bool CalcEstatics::EvaluateFrame(QMTopology *top) {
    cout << "Doing Estatics\n";
    vector<CrgUnit *> lcharges = top->CrgUnits();
    Topology atop;
    atop.setBox(top->getBox());
    for (vector<CrgUnit *>::iterator itl = lcharges.begin(); itl != lcharges.end(); itl++) {
        top->AddAtomisticBeads(*itl, &atop);
    }
    cout << "Number of charge units in top " << lcharges.size() << endl;
    cout << "Number of molecules in atop " << atop.MoleculeCount() << endl;
    MoleculeContainer::iterator imol;
    for (imol = atop.Molecules().begin(); imol != atop.Molecules().end(); imol++) {
        Molecule *mol = *imol;

        //Compute the EstaticPotentialEnergy for molecule mol if it is neutral
        double neutr = CalcPot(&atop, mol);

        CrgUnit *crg = mol->getUserData<CrgUnit > ();

        //Change the charges on the molecule mol to the occupied=charged state
        top->CopyChargesOccupied(crg, mol);

        //Compute the EstaticPotentialEnergy for molecule mol if it is charged
        double crged = CalcPot(&atop, mol);

        //If the site energy is large and positive (crged >> neutr) molecule mol wants to stay neutral
        //This is consistent with marcus_rates.h since w_12 (from 1 to 2) contains (dG_12 + \lambda)^2 and dG_12=e_2-e_1
        crg->setEnergy(crged - neutr);

        //Copy back the charges as if molecule mol would be neutral so that the loop can pick up the next molecule
        top->CopyCharges(crg, mol);
        cout << "Estatic energy [eV] for charged / neutral / crg-neutr=espilon: " << crged << " " << neutr << " " << crged - neutr << "\n";
    }

    atop.Cleanup();

    return true;
}

//Calculate Estatics

double CalcEstatics::CalcPot(Topology *atop, Molecule *mol) //wegen Übergabe per * unten ->
{
    //relative dielectric constant
    double epsilon_dielectric = 3.5;
    //estatic energy including contributions from all other molecules in eV
    double pot = 0.0;
    //Coulombs constant including conversion factors (elementary charge and nm) so that pot is in eV
    double k_c = 1.602176487e-19 / (1.0e-9 * 8.854187817e-12 * 4.0 * M_PI);
    cout << "Doing Estatics with epsilon_relative=" << epsilon_dielectric << endl;

    //Do loop over all other molecules imol that are not mol
    MoleculeContainer::iterator imol;
    for (imol = atop->Molecules().begin(); imol != atop->Molecules().end(); imol++) {
        if (*imol == mol) continue;
        //Check whether PBC have to be taken into account
        vec s;
        //We should replace getUserData by a funtion GetCom for molecules
        vec bcs = atop->BCShortestConnection(mol->getUserData<CrgUnit > ()->GetCom(), (*imol)->getUserData<CrgUnit > ()->GetCom());
        vec dist = (*imol)->getUserData<CrgUnit > ()->GetCom() - mol->getUserData<CrgUnit > ()->GetCom();
        vec diff = bcs - dist;
        
        s = diff;
        if (abs(dist) > abs(bcs)) {
            s.setX(diff.x());
            s.setY(diff.y());
            s.setZ(diff.z());
        } else {
            s.setX(0.0);
            s.setY(0.0);
            s.setZ(0.0);        //if abs(dist)< abs(bcs) something is wrong. This happens in triclinic boxes.
            cout << dist << "; " << bcs << abs(dist) << ";" << abs(bcs) << endl;
        }

        //cout<<" s="<<s<<endl;
        for (int i = 0; i < mol->BeadCount(); i++)
            for (int j = 0; j != (*imol)->BeadCount(); j++) {
                Bead *bi = mol->getBead(i);
                Bead *bj = (*imol)->getBead(j);
                //distance vector may have to be shifted by s
                vec r_v = bi->getPos()-(bj->getPos() + s);
                double r = abs(r_v);
                double qi = bi->getQ();
                double qj = bj->getQ();
                //cout<<bj->getPos()<<r_v<<r<<qi<<qj<<epsilon_dielectric<<pot<<qi*qj/epsilon_dielectric*1/r<<endl;
                pot += k_c * qi * qj / (epsilon_dielectric * r); //in eV
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
