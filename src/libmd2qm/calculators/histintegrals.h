/*
 * File:   calc_integrals.h
 * Author: schrader
 *
 * Created on July 7, 2010
 */

#ifndef _HIST_INTEGRALS_H
#define	_HIST_INTEGRALS_H

#include "qmpair.h"
#include "paircalculator.h"
#include <votca/tools/histogramnew.h>

class CalcHistIntegrals : public PairCalculator
{
public:
    CalcHistIntegrals() {};
    ~CalcHistIntegrals() {};

    void Initialize(QMTopology *top, Property *options);
    void EvaluatePair(QMTopology *top, QMPair *pair);
    void EndEvaluate(QMTopology *top);

private:
    HistogramNew histogram;
    double _min, _max;
    int _nbins;
    string _outfile;
};

inline void CalcHistIntegrals::Initialize(QMTopology *top, Property *options){
    _min = options->get("options.histintegrals.min").as<double>();
    _max = options->get("options.histintegrals.max").as<double>();
    _nbins = options->get("options.histintegrals.nbins").as<int>();
    _outfile = options->get("options.histintegrals.file").as<string>();

    histogram.Initialize(_min,_max,_nbins);
}

inline void CalcHistIntegrals::EndEvaluate(QMTopology *top){
    cout << "Writing Distribution of Log10(J_eff^2) to file " << _outfile << endl;
    histogram.data().Save(_outfile);
}

inline void CalcHistIntegrals::EvaluatePair(QMTopology *top, QMPair *pair){
    double logJ2eff = log10(pair->calcJeff2());
    histogram.Process( logJ2eff );
}

#endif	/* _HIST_INTEGRALS_H */