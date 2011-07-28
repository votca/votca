#include <stdlib.h>
#include <math.h>
#include <list>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include "ecoulomb.h"

void Ecoulomb::Initialize(QMTopology *top, Property *options) {
    _options = options;
    if (options->exists("options.ecoulomb.method")) {
        if (options->get("options.ecoulomb.method").as<string > () == "distdep") {
            if (options->exists("options.ecoulomb.epsilon")) {
                _epsilon_dielectric = options->get("options.ecoulomb.epsilon").as<double>();
            } else {
                //_epsilon_dielectric = 3.0;
                std::runtime_error("Error in ecoulomb: dielectric constant is not provided");
            }
            if (options->exists("options.ecoulomb.screening")) {
                _s_eps = options->get("options.ecoulomb.screening").as<double>();
            } else {
                //_s_eps = 3.0;
                std::runtime_error("Error in ecoulomb: screening length is not provided.");
            }
            cout << "Doing distance-dependent screening with eps " << _epsilon_dielectric << " and screening " << _s_eps << endl;
	    _estatic_method = &Ecoulomb::dist_dep_eps;
        }
        else if (options->get("options.ecoulomb.method").as<string > () == "simple") {
             if (options->exists("options.ecoulomb.epsilon")) {
                 _epsilon_dielectric = options->get("options.ecoulomb.epsilon").as<double>();
             } else {
                  //_epsilon_dielectric = 3.0;
                  std::runtime_error("Error in ecoulomb: dielectric constant is not provided.");
             }
             cout << "Doing simple estatic with constant epsilon = " << _epsilon_dielectric << endl;
	     _estatic_method = &Ecoulomb::constant_epsilon;
        }

        else throw std::runtime_error("Error in ecoulomb: the method (simple or distdep) is not specified.");
    } else throw std::runtime_error("Error in ecoulomb: the method (simple or distdep) is not specified.");
}

bool Ecoulomb::EvaluateFrame(QMTopology *top) {
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
        QMCrgUnit *crg = mol->getUserData<QMCrgUnit > ();
        //Compute the Coulomb sum of the system for a neutral molecule mol
        neutr = CalcPot(&atop, mol);
        //Change the partial charges of mol to the charged state
        top->CopyChargesOccupied(crg, mol);
        //Compute the Coulommb sum of the system for mol in a charged state
        crged = CalcPot(&atop, mol);
        //Copy back the charges as if molecule mol would be neutral so that the loop can pick up the next molecule
        top->CopyCharges(crg, mol);
	//If the site energy is large and positive (crged >> neutr) molecule mol wants to stay neutral
	//This is consistent with marcus_rates.h since w_12 (from 1 to 2) contains (dG_12 + \lambda)^2 and dG_12=e_2-e_1
        //Attention not to overwrite the site energies from list_charges.xml
        crg->setDouble("energy_coulomb",crged - neutr);
        //cout << "Estatic energy [eV] for charged / neutral / crg-neutr=espilon: " << crged << " " << neutr << " " << crged - neutr << "\n";
        cout << "Ecoulomb for crgunit " << crg->getId() << " at pos " << crg->GetCom() << " is " << crged - neutr << " eV\n";
    }   
    atop.Cleanup();
    return true;
}

// Couloumb interactions (sum over all pairs) 
double Ecoulomb::CalcPot(Topology *atop, Molecule *mol) 
{
    //estatic energy including contributions from all other molecules in eV
    double pot = 0.0;
    //Coulomb constant including conversion factors (elementary charge and nm). Resulting potential is in eV
    double k_c = 1.602176487e-19 / (1.0e-9 * 8.854187817e-12 * 4.0 * M_PI);
    //Do loop over all other molecules imol that are not mol
    MoleculeContainer::iterator imol;
    for (imol = atop->Molecules().begin(); imol != atop->Molecules().end(); imol++) {
        if (*imol == mol) continue;
        // periodic boundary-corrected distance between the ceneters of mass
	// used for the distance-dependent screening (we should replace getUserData by a funtion GetCom for molecules)
	vec bcs = atop->BCShortestConnection(mol->getUserData<CrgUnit > ()->GetCom(), (*imol)->getUserData<CrgUnit > ()->GetCom());
        vec dist = (*imol)->getUserData<CrgUnit > ()->GetCom() - mol->getUserData<CrgUnit > ()->GetCom();
        vec diff = bcs - dist;

	for (int i = 0; i < mol->BeadCount(); i++)
            for (int j = 0; j != (*imol)->BeadCount(); j++) {
                Bead *bi = mol->getBead(i);
                Bead *bj = (*imol)->getBead(j);
                // distance between two charges taking into account periodic boundary conditions
                vec r_v = bi->getPos()-(bj->getPos() + diff);
                pot += k_c * bi->getQ() * bj->getQ() / ( (this->*_estatic_method)(abs(bcs)) * abs(r_v) ); //in eV
            }
    }
    return pot;

}

double Ecoulomb::dist_dep_eps(const double &dist) {
    return _epsilon_dielectric - (_epsilon_dielectric - 1.0) * exp(-_s_eps * dist)*(1.0 + _s_eps * dist + 0.5 * _s_eps * dist * _s_eps * dist);    
}

double Ecoulomb::constant_epsilon(const double &dist) {
    return _epsilon_dielectric;
}
