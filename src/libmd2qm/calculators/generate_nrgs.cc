/* 
 * File:   generate_nrgs.cc
 * Author: lukyanov
 * 
 * Created on July 21, 2010, 4:46 PM
 */

#include <stdlib.h>
#include "generate_nrgs.h"
#include <math.h>
#include <list>
#include <votca/tools/random.h>
#include "qmnblist.h"
#include <votca/tools/average.h>
#include "votca/csg/nblistgrid.h"

void GenerateNrgs::Initialize(QMTopology *top, Property *options) {

   // read in _sigma
   if (options->exists("options.generate_energies.sigma")) {
   	_sigma = options->get("options.generate_energies.sigma").as<double>();
   }
   else {
	_sigma = 1.0;
	cout << "Warning: sigma of site energy distribution is not provided, using default "
             << "sigma = " << _sigma << endl;
   }

   // read in _correl
   if (options->exists("options.generate_energies.correl")) {
   	_correl = options->get("options.generate_energies.correl").as<bool>();
   }
   else {
	_correl = false;
	cout << "Warning: correlations of site energies are not specified, using non-correlated" << endl;
   }

   // print info
   if (_correl) {
            // read in _cutoff
            if (options->exists("options.generate_energies.cutoff")) {
                _cutoff = options->get("options.generate_energies.cutoff").as<double>();
            }
            else {
                _cutoff = 1.5;
                cout << "Warning: cutoff for correlations is not specified, using default" << endl;
            }
        cout << "Generating correlated site energies with sigma = " << _sigma << " and cutoff = " << _cutoff << endl;
   }
   else {
        cout << "Generating non-correlated site energies with sigma = " << _sigma << endl;
   }

    /// Initialize the random number generator
    Random::init(14, 122, 472, 1912);
}




bool GenerateNrgs::EvaluateFrame(QMTopology *top) {

    if (_correl) {
        // clear _tmp_energy to avoid accumulation of stuff
        // from the previous frame
        _tmp_energy.clear();
        AssignCorrelated( top );
    }
    else {   
        AssignGaussian( top );
    }

    return true;
}

void GenerateNrgs::AssignGaussian(QMTopology *top) {

    vector<CrgUnit *> lcharges = top->CrgUnits();
    vector<CrgUnit *>::iterator itl;

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        (*itl)->setEnergy( Random::rand_gaussian(_sigma) );
    }

}

void GenerateNrgs::AssignCorrelated(QMTopology *top) {
    // First assign gaussian energies
    AssignGaussian( top );

    // working topology = topology to do work
    Topology mytop;
    BeadList list1;
    NBListGrid mynbl;

    mytop.setBox(top->getBox());

    vector<CrgUnit *> lcharges = top->CrgUnits();
    vector<CrgUnit *>::iterator itl;

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        // make a simple spherical bead for a new topology
        BeadType *tmptype = mytop.GetOrCreateBeadType("no");
        Bead * tmpbead = mytop.CreateBead(1, "bead", tmptype, 0, 1, 0);
        // position of the bead = position of the corresponding crgunit
        tmpbead->setPos( (*itl)->GetCom() );
        // additionally it stores a pointer to the corresponding crgunit
        tmpbead->setUserData( (*itl) );
        
    }

    list1.Generate(mytop, "*");
    mynbl.SetMatchFunction(this, &GenerateNrgs::MyMatchingFunction);
    mynbl.setCutoff(_cutoff);
    mynbl.Generate( list1, false );

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        // e1+e2+...+eN/ sqrt(N) looks weird, but should be done to get the same sigma
        (*itl)->setEnergy( _tmp_energy[(*itl)].getAvg() * sqrt( _tmp_energy[(*itl)].getN() ) );
    }

}

bool GenerateNrgs::MyMatchingFunction(Bead *bead1, Bead *bead2, const vec & r) {

    CrgUnit *crg1 = bead1->getUserData<CrgUnit>();
    CrgUnit *crg2 = bead2->getUserData<CrgUnit>();


    double e1 = crg1->getEnergy();
    double e2 = crg2->getEnergy();

    _tmp_energy[crg1].Process(e2);
    _tmp_energy[crg2].Process(e1);
    
    return false;
}

/*
     QMNBList  &nblist = top->nblist();
    vector<CrgUnit *> lcharges = top->CrgUnits();
    std::map<CrgUnit *, Average<double> > tmp_energy;

    vector<CrgUnit *>::iterator itl;
    QMNBList::iterator piter;

    // loop over all pairs
    for (piter = nblist.begin(); piter!=nblist.end(); ++piter) {
        double e1 = (*piter)->first->getEnergy();
        double e2 =  (*piter)->second->getEnergy();

        tmp_energy[(*piter)->first].Process(e2);
        tmp_energy[(*piter)->second].Process(e1);
    }

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        // e1+e2+...+eN/ sqrt(N) looks weird, but should be done to get the same sigma
        (*itl)->setEnergy( tmp_energy[(*itl)].getAvg() * sqrt( tmp_energy[(*itl)].getN() ) );
    }
 
 */