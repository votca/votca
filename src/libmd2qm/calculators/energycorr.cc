/* 
 * File:   energycorr.cc
 * Author: lukyanov
 * 
 * Created on July 28, 2010, 6:03 PM
 */

#include <stdlib.h>
#include "energycorr.h"
#include <math.h>
#include <list>
#include "votca/tools/average.h"


void EnergyCorr::Initialize(QMTopology* top, Property* options) {

   // read in _nbins
   if (options->exists("options.energy_corr.nbins")) {
   	_nbins = options->get("options.energy_corr.nbins").as<int>();
   }
   else {
	_nbins = 15;
	cout << "Warning: nbins for site energy distribution is not provided, using default "
             << "nbins = " << _nbins << endl;
   }

   // read in min and max
   _min = options->get("options.energy_corr.min").as<double>();
   _max = options->get("options.energy_corr.max").as<double>();

   // read in _outfile
   if (options->exists("options.energy_corr.file")) {
   	_outfile = options->get("options.energy_corr.file").as<string>();
   }
   else {
	_outfile = "energy_corr.dat";
	cout << "Warning: output filename for site energy distribution is not provided, using default "
             << "filename = " << _outfile << endl;
   }

   // initialize histograms
   _hist_corr.Initialize(_min, _max, _nbins);
   _hist_count.Initialize(_min, _max, _nbins);
   _hist_corr2.Initialize(_min, _max, _nbins);

   // set yerror = true
   _hist_corr.data().SetHasYErr(true);
}

bool EnergyCorr::EvaluateFrame(QMTopology* top) {

    double mean_energy, stdev_energy;
    Average<double> energy_av;

    vector<CrgUnit *> lcharges = top->CrgUnits();
    vector<CrgUnit *>::iterator itl;

    QMNBList  &nblist = top->nblist();
    QMNBList::iterator piter;

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        energy_av.Process ( (*itl)->getEnergy() );
    }

    // mean energy and standard deviation for a particular snapshot
    mean_energy = energy_av.getAvg();
    stdev_energy = energy_av.CalcSig2();

    // loop over all pairs
    for (piter = nblist.begin(); piter!=nblist.end(); ++piter) {
        double e1 = (*piter)->first->getEnergy();
        double e2 =  (*piter)->second->getEnergy();
        double dist = (*piter)->dist();

        double contrib = (e1-mean_energy)*(e2-mean_energy)/stdev_energy;
        
        _hist_count.Process( dist );
        _hist_corr.Process( dist, contrib );
        _hist_corr2.Process( dist, contrib*contrib );
    }
    
    return true;
}

void EnergyCorr::EndEvaluate(QMTopology* top) {
    // normalize
    _hist_corr.data().y() = element_div( _hist_corr.data().y(), _hist_count.data().y() );
    _hist_corr2.data().y() = element_div( _hist_corr2.data().y(), _hist_count.data().y() );

    // calculate error bars: UGLY!!!
    for (int i=0; i<_hist_corr.data().yerr().size(); i++) {
        _hist_corr.data().yerr(i) = sqrt ( (_hist_corr2.data().y(i) -
                _hist_corr.data().y(i)*_hist_corr.data().y(i) ) / _hist_count.data().y(i) );
    }
    
    // Save to file
    _hist_corr.data().Save(_outfile);
}
