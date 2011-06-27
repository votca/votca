#ifndef _HIST_ENERGETICDISORDER_H
#define	_HIST_ENERGETICDISORDER_H

#include "qmpair.h"
#include "paircalculator.h"
#include <votca/tools/histogramnew.h>

/**
	\brief Compute a histogram of site energy differences from the neighborlist

Callname: histenergeticdisorder
*/
class CalcHistEnergeticDisorder : public PairCalculator
{
public:
    CalcHistEnergeticDisorder() {};
    ~CalcHistEnergeticDisorder() {};

    const char *Description() { return "Compute a histogram of site energy differences from the neighborlist"; }

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
    cout << "Writing Distribution of site energy differences to file " << _outfile << endl;
    histogram.data().Save(_outfile);
}

inline void CalcHistEnergeticDisorder::EvaluatePair(QMTopology *top, QMPair *pair){
    double energy1 = pair->Crg1()->getEnergy();
    double energy2 = pair->Crg2()->getEnergy();
    histogram.Process( energy1-energy2 );
}

#endif	/* _HIST_ENERGETICDISORDER_H */
