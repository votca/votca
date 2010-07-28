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
        cout << "Generating correlated site energies with sigma = " << _sigma << endl;
   }
   else {
        cout << "Generating non-correlated site energies with sigma = " << _sigma << endl;
   }

    /// Initialize the random number generator
    Random::init(14, 122, 472, 1912);
}




bool GenerateNrgs::EvaluateFrame(QMTopology *top) {

    if (_correl) {
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
}

