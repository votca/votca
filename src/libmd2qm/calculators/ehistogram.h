#ifndef _EHISTOGRAM_H
#define	_EHISTOGRAM_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/paircalculator.h>
#include <votca/tools/histogramnew.h>

/**
	\brief Histogram of site energy differences of neighbor list pairs

Callname: ehistogram
*/
class Ehistogram : public PairCalculator
{
public:
    Ehistogram() {};
    ~Ehistogram() {};

    const char *Description() { return "Histogram of site energy differences of neighbor list pairs"; }

    void Initialize(QMTopology *top, Property *options);
    void EvaluatePair(QMTopology *top, QMPair *pair);
    void EndEvaluate(QMTopology *top);

private:
    HistogramNew histogram;
    double _min, _max;
    int _nbins;
    string _outfile;
};

inline void Ehistogram::Initialize(QMTopology *top, Property *options){
    _min = options->get("options.ehistogram.min").as<double>();
    _max = options->get("options.ehistogram.max").as<double>();
    _nbins = options->get("options.ehistogram.nbins").as<int>();
    _outfile = options->get("options.ehistogram.file").as<string>();

    histogram.Initialize(_min,_max,_nbins);
}

inline void Ehistogram::EndEvaluate(QMTopology *top){
    cout << "Writing Distribution of site energy differences to file " << _outfile << endl;
    histogram.data().Save(_outfile);
}

inline void Ehistogram::EvaluatePair(QMTopology *top, QMPair *pair){
    double energy1 = pair->first->getTotalEnergy();
    double energy2 = pair->second->getTotalEnergy();
    histogram.Process( energy1-energy2 );
}

#endif	/* _EHISTOGRAM_H */
