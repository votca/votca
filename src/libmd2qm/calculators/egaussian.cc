#include <stdlib.h>
#include "egaussian.h"
#include <math.h>
#include <list>
#include <votca/tools/random.h>
#include <votca/ctp/qmnblist.h>
#include <votca/tools/average.h>
#include <votca/csg/nblistgrid.h>

void Egaussian::Initialize(QMTopology *top, Property *options) {

   // read in _sigma
   if (options->exists("options.egaussian.sigma")) {
   	_sigma = options->get("options.egaussian.sigma").as<double>();
   }
   else {
	//_sigma = 1.0;
	std::runtime_error("Error in egaussian: variance (sigma) of site energy distribution is not provided");
   }

   // read in _correl
   if (options->exists("options.egaussian.correlation")) {
   	_correl = options->get("options.generate_energies.correl").as<bool>();
   }
   else {
	_correl = false;
	cout << "Warning: correlation of site energies is not specified, using non-correlated" << endl;
   }

   // print info
   if (_correl) {
            // read in _cutoff
            if (options->exists("options.egaussian.cutoff")) {
                _cutoff = options->get("options.generate_energies.cutoff").as<double>();
            }
            else {
                //_cutoff = 1.5;
                std::runtime_error("Error in egaussian: cutoff for correlations is not specified, using default");
            }
        cout << "Generating correlated site energies with sigma = " << _sigma << " and cutoff = " << _cutoff << endl;
   }
   else {
        cout << "Generating non-correlated site energies with sigma = " << _sigma << endl;
   }

    /// Initialize the random number generator
    Random::init(14, 122, 472, 1912);
}




bool Egaussian::EvaluateFrame(QMTopology *top) {

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

void Egaussian::AssignGaussian(QMTopology *top) {

    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    vector<QMCrgUnit *>::iterator itl;

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        (*itl)->setDouble("energy_coulomb", Random::rand_gaussian(_sigma) );
    }
}

void Egaussian::AssignCorrelated(QMTopology *top) {
    // First assign gaussian energies
    AssignGaussian( top );

    // working topology = topology to do work
    Topology mytop;
    BeadList list1;
    NBListGrid mynbl;

    mytop.setBox(top->getBox());

    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    vector<QMCrgUnit *>::iterator itl;

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
    mynbl.SetMatchFunction(this, &Egaussian::MyMatchingFunction);
    mynbl.setCutoff(_cutoff);
    mynbl.Generate( list1, false );

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        // e1+e2+...+eN/ sqrt(N) - normalization to get the same sigma
        (*itl)->setDouble("energy_coulomb", _tmp_energy[(*itl)].getAvg() * sqrt( _tmp_energy[(*itl)].getN() ) );
    }

}

bool Egaussian::MyMatchingFunction(Bead *bead1, Bead *bead2, const vec & r, const double notused) {

    QMCrgUnit *crg1 = bead1->getUserData<QMCrgUnit>();
    QMCrgUnit *crg2 = bead2->getUserData<QMCrgUnit>();


    double e1 = crg1->getTotalEnergy();
    double e2 = crg2->getTotalEnergy();

    _tmp_energy[crg1].Process(e2);
    _tmp_energy[crg2].Process(e1);
    
    return false;
}
