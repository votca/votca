/*
 * File:   calc_estatics.cpp
 * Author: mayfalk
 *
 * Created on May 21, 2010, 11:21 PM
 */


#include <stdlib.h>
#include "estatics.h"
#include <math.h>
#include <list>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>

void CalcEstatics::Initialize(QMTopology *top, Property *options) {
    _options = options;
    if (options->exists("options.estatic_params.epsilon_dielectric")) {
        _epsilon_dielectric = options->get("options.estatic_params.epsilon_dielectric").as<double>();
    } else {
        _epsilon_dielectric = 3.0;
        cout << "Warning: dielectric constant in not provided, using default" << endl;
    }

    if (options->exists("options.estatic_params.s_eps")) {
        _s_eps = options->get("options.estatic_params.s_eps").as<double>();
    } else {
        _s_eps = 3.0;
        cout << "Warning: s_eps is not provided, using default" << endl;
    }

    if (options->exists("options.estatic_params.estatic_method")) {
        if (options->get("options.estatic_params.estatic_method").as<string > () == "distdep") {
            _estatic_method = &CalcEstatics::dist_dep_eps;
            cout << "Doing distance-dependent-eps estatic with eps " << _epsilon_dielectric << " and s_eps " << _s_eps << endl;
        }
        if (options->get("options.estatic_params.estatic_method").as<string > () == "simple") {
            _estatic_method = &CalcEstatics::constant_epsilon;
            cout << "Doing simple estatic with eps " << _epsilon_dielectric << endl;
        } else if (options->get("options.estatic_params.estatic_method").as<string > () == "dipoles") {
            _estatic_method = &CalcEstatics::constant_epsilon;
            cout << "Doing interactions including dipoles with eps =1 " << endl;
        } else throw std::runtime_error("Error in CalcEstatics::Initialize : no such estatic method, should be simple or distdep");
    } else throw std::runtime_error("Error in CalcEstatics:Initialize : no estatic_method specified");
}

bool CalcEstatics::EvaluateFrame(QMTopology *top) {
    vector<CrgUnit *> lcharges = top->CrgUnits();
    Topology atop;
    atop.setBox(top->getBox());
    for (vector<CrgUnit *>::iterator itl = lcharges.begin(); itl != lcharges.end(); itl++) {
        top->AddAtomisticBeads(*itl, &atop);
    }
    cout << "Number of charge units in top " << lcharges.size() << endl;
    cout << "Number of molecules in atop " << atop.MoleculeCount() << endl;
    double crged, neutr;

    MoleculeContainer::iterator imol;
    for (imol = atop.Molecules().begin(); imol != atop.Molecules().end(); imol++) {
        Molecule *mol = *imol;
        CrgUnit *crg = mol->getUserData<CrgUnit > ();
        //if induced dipoles then directly go to charged system
        if (_options->get("options.estatic_params.estatic_method").as<string > () == "dipoles") {

            //Change the charges on the molecule mol to the occupied=charged state
            top->CopyChargesOccupied(crg, mol);
            //read in the induced dipole moments
            string s;
            int nr = mol->getId();
            string filename = "d_" + boost::lexical_cast<string > (nr);
            ifstream indata;
            indata.open(filename.c_str());
            for (BeadContainer::iterator iter = atop.Beads().begin();
                    iter != atop.Beads().end(); ++iter) { // keep reading until end-of-file
                Bead *b = *iter;
                getline(indata, s);
                Tokenizer tok(s, " ");
                vector<double> v;
                tok.ConvertToVector(v);
                if (v.size() != 3)
                    throw std::ios_base::failure("invalid dipole format. should be dx dy dz in Debeye separated by one space");
                double vx = v[0];
                double vy = v[1];
                double vz = v[2];
                b->setUserData<vec > (new vec(vx, vy, vz));
            }
            indata.close();
            crged = CalcPot_Dipole(&atop, mol);
            neutr = 0.0;
            top->CopyCharges(crg, mol);
        }
            //otherwise do the regular estatics shit
        else {
            //Compute the EstaticPotentialEnergy for molecule mol if it is neutral
            neutr = CalcPot(&atop, mol);

            CrgUnit *crg = mol->getUserData<CrgUnit > ();

            //Change the charges on the molecule mol to the occupied=charged state
            top->CopyChargesOccupied(crg, mol);

            //Compute the EstaticPotentialEnergy for molecule mol if it is charged
            crged = CalcPot(&atop, mol);
            //Copy back the charges as if molecule mol would be neutral so that the loop can pick up the next molecule
            top->CopyCharges(crg, mol);
        }

        //If the site energy is large and positive (crged >> neutr) molecule mol wants to stay neutral
        //This is consistent with marcus_rates.h since w_12 (from 1 to 2) contains (dG_12 + \lambda)^2 and dG_12=e_2-e_1
        crg->setEnergy(crged - neutr);

        //        cout << "Estatic energy [eV] for charged / neutral / crg-neutr=espilon: " << crged << " " << neutr << " " << crged - neutr << "\n";
        cout << "estatic energy [eV] for crgunit " << crg->getId() << " at pos " << crg->GetCom() << " is " << crged - neutr << "\n";

        if (_options->get("options.estatic_params.estatic_method").as<string > () == "dipoles") {
            for (BeadContainer::iterator iter = atop.Beads().begin();
                    iter != atop.Beads().end(); ++iter) {
                Bead *b = *iter;
               // cout << "if mol " << crg->getId() << "at pos " << crg->GetCom() << " is charged, Bead " << b->getId() << " at pos " << b->getPos() << " has charge " << b->getQ() << " and dipole xyz [D]" << *b->getUserData<vec > () << "\n";
                delete b->getUserData<vec > ();
                b->setUserData<void>(NULL);
            }
        }
    }

    atop.Cleanup();

    return true;
}

//Calculate Estatics

double CalcEstatics::CalcPot(Topology *atop, Molecule *mol) //wegen Übergabe per * unten ->
{
    //estatic energy including contributions from all other molecules in eV
    double pot = 0.0;
    //Coulombs constant including conversion factors (elementary charge and nm) so that pot is in eV
    double k_c = 1.602176487e-19 / (1.0e-9 * 8.854187817e-12 * 4.0 * M_PI);

    //Do loop over all other molecules imol that are not mol
    MoleculeContainer::iterator imol;
    for (imol = atop->Molecules().begin(); imol != atop->Molecules().end(); imol++) {
        if (*imol == mol) continue;
        //Check whether PBC have to be taken into account
        //We should replace getUserData by a funtion GetCom for molecules
        vec bcs = atop->BCShortestConnection(mol->getUserData<CrgUnit > ()->GetCom(), (*imol)->getUserData<CrgUnit > ()->GetCom());
        vec dist = (*imol)->getUserData<CrgUnit > ()->GetCom() - mol->getUserData<CrgUnit > ()->GetCom();
        vec diff = bcs - dist;

        //cout<<" s="<<s<<endl;
        for (int i = 0; i < mol->BeadCount(); i++)
            for (int j = 0; j != (*imol)->BeadCount(); j++) {
                Bead *bi = mol->getBead(i);
                Bead *bj = (*imol)->getBead(j);
                //distance vector may have to be shifted by s
                vec r_v = bi->getPos()-(bj->getPos() + diff);
                double r = abs(r_v);
                double qi = bi->getQ();
                double qj = bj->getQ();
                //cout<<bj->getPos()<<r_v<<r<<qi<<qj<<_epsilon_dielectric<<pot<<qi*qj/_epsilon_dielectric*1/r<<endl;
                //pot += k_c * qi * qj / ( (this->*_estatic_method)(r) * r ); //in eV
                pot += k_c * qi * qj / ((this->*_estatic_method)(abs(bcs)) * r); //in eV
            }
    }
    return pot;

}

double CalcEstatics::CalcPot_Dipole(Topology *atop, Molecule *mol) //wegen Übergabe per * unten ->
{
    //estatic energy including contributions from all other molecules in eV
    double pot;
    double pot_qq=0.0;
    double pot_qd=0.0;
    double pot_dd=0.0;
    //Coulombs constant including conversion factors (elementary charge and nm) so that pot is in eV
    double k_cc = 1.602176487e-19 / (1.0e-9 * 8.854187817e-12 * 4.0 * M_PI);
    double k_cd = k_cc * 0.02082; //1 Debeye is 3.33564 e-30 C*m or 0.02082*e*nm
    double k_dd = k_cc * 0.02082*0.02082;
    //Now we have to include interactions of the other molecules among themselves
    //Do loop over all molecules cmol
    MoleculeContainer::iterator cmol;
    //Do loop over all other molecules imol that are not cmol
    MoleculeContainer::iterator imol;
    for (cmol = atop->Molecules().begin(); cmol != atop->Molecules().end(); cmol++) {
        for (imol = atop->Molecules().begin(); imol != atop->Molecules().end(); imol++) {
            //only count interaction if mol_id of c is smaller than mol_id of i
            if ((*imol)->getId() <= (*cmol)->getId()) continue;
            //Check whether PBC have to be taken into account
            //We should replace getUserData by a funtion GetCom for molecules
            vec bcs = atop->BCShortestConnection((*cmol)->getUserData<CrgUnit > ()->GetCom(), (*imol)->getUserData<CrgUnit > ()->GetCom());
            vec dist = (*imol)->getUserData<CrgUnit > ()->GetCom() - (*cmol)->getUserData<CrgUnit > ()->GetCom();
            vec diff = bcs - dist;

            //cout<<" s="<<s<<endl;
            for (int i = 0; i != (*cmol)->BeadCount(); i++)
                for (int j = 0; j != (*imol)->BeadCount(); j++) {
                    Bead *bi = (*cmol)->getBead(i);
                    Bead *bj = (*imol)->getBead(j);
                    //distance vector may have to be shifted by s
                    vec r_v = bi->getPos()-(bj->getPos() + diff); //ri-rj so that r_v points from j to i
                    double r = abs(r_v);
                    double r3 = r * r*r;
                    double qi = bi->getQ();
                    double qj = bj->getQ();
                    vec di = *bi->getUserData<vec > ();
                    vec dj = *bj->getUserData<vec > ();
                    //charge charge interaction in eV
                    pot_qq += k_cc * qi * qj / r;
                    //dipol-charge interaction in eV between charge on i and dipole on j
                    pot_qd += k_cd * qi * (dj.getX() * r_v.getX() + dj.getY() * r_v.getY() + dj.getZ() * r_v.getZ()) / r3;
                    //dipol-charge interaction in eV if charge on j, dipole on i. Note that r_v points from j to i, so -1.0!
                    pot_dd += k_dd * (-1.0) * qj * (di.getX() * r_v.getX() + di.getY() * r_v.getY() + di.getZ() * r_v.getZ()) / r3;
                    //dipole-dipole  interaction
                    pot_dd += k_dd / r3 * (dj.getX() * di.getX() + dj.getY() * di.getY() + dj.getZ() * di.getZ() + 3.0 * (di.getX() * r_v.getX() + di.getY() * r_v.getY() + di.getZ() * r_v.getZ())*(dj.getX() * r_v.getX() + dj.getY() * r_v.getY() + dj.getZ() * r_v.getZ()));
                    
                }
        }
    }
    pot=pot_qq+pot_qd+pot_dd;
    cout << "Q-Q: "<<pot_qq << "    Q-D: "<<pot_qd<< "  D-D: "<< pot_dd<<"  total: "<< pot <<"\n";
    return pot;
}

double CalcEstatics::dist_dep_eps(const double &dist) {
    return _epsilon_dielectric - (_epsilon_dielectric - 1.0) * exp(-_s_eps * dist)*(1.0 + _s_eps * dist + 0.5 * _s_eps * dist * _s_eps * dist);
}

double CalcEstatics::constant_epsilon(const double &dist) {
    return _epsilon_dielectric;
}
