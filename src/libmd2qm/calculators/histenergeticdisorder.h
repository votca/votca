/*
 * File:   calc_integrals.h
 * Author: schrader
 *
 * Created on March 28, 2011
 */

#ifndef _HIST_ENERGETICDISORDER_H
#define	_HIST_ENERGETICDISORDER_H

#include "qmpair.h"
#include "paircalculator.h"
#include <votca/tools/histogramnew.h>

class CalcHistEnergeticDisorder : public PairCalculator
{
public:
    CalcHistEnergeticDisorder() {};
    ~CalcHistEnergeticDisorder() {};

    void Initialize(QMTopology *top, Property *options);
    void EvaluatePair(QMTopology *top, QMPair *pair);
    void EndEvaluate(QMTopology *top);

private:
    HistogramNew histogram;
    double _min, _max;
    int _nbins;
    string _outfile;
};

inline void CalcHistEnergeticDisorder::Initialize(QMTopology *top, Property *options){
    _min = options->get("options.histenergeticdisorder.min").as<double>();
    _max = options->get("options.histenergeticdisorder.max").as<double>();
    _nbins = options->get("options.histenergeticdisorder.nbins").as<int>();
    _outfile = options->get("options.histenergeticdisorder.file").as<string>();

    histogram.Initialize(_min,_max,_nbins);
}

inline void CalcHistEnergeticDisorder::EndEvaluate(QMTopology *top){
    cout << "Writing Distribution of energetic disorder to file " << _outfile << endl;
    histogram.data().Save(_outfile);
}

inline void CalcHistEnergeticDisorder::EvaluatePair(QMTopology *top, QMPair *pair){
    double energy1 = pair->Crg1()->getEnergy();
    double energy2 = pair->Crg2()->getEnergy();
    histogram.Process( energy1-energy2 );
}

#endif	/* _HIST_ENERGETICDISORDER_H */
