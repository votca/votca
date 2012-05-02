#ifndef IANALYZE_H
#define IANALYZE_H

#include <votca/ctp/qmcalculator.h>
#include <math.h>


namespace votca { namespace ctp {

class IAnalyze : public QMCalculator
{
public:

    string  Identify() { return "IAnalyze"; }

    void    Initialize(Topology *top, Property *options);
    bool    EvaluateFrame(Topology *top);
    void    IHist(Topology *top, int state);

private:

    double      _resolution_logJ2;
    vector<int> _states;

};


void IAnalyze::Initialize(Topology *top, Property *opt) {

    string key = "options.ianalyze";

    _resolution_logJ2 = opt->get(key+".resolution_logJ2").as< double >();
    _states = opt->get(key+".states").as< vector<int> >();
}


bool IAnalyze::EvaluateFrame(Topology *top) {

    QMNBList &nblist = top->NBList();

    if (!nblist.size()) {
        cout << endl << "... ... No pairs in topology. Skip...";
        return 0;
    }

    for (int i = 0; i < _states.size(); ++i) {
        this->IHist(top, _states[i]);
    }
}


void IAnalyze::IHist(Topology *top, int state) {

    QMNBList &nblist = top->NBList();
    QMNBList::iterator nit;

    double MIN = log10(nblist.front()->getJeff2(state));
    double MAX = log10(nblist.front()->getJeff2(state));

    // Collect J2s from pairs
    vector< double > J2s;
    J2s.reserve(nblist.size());

    for (nit = nblist.begin(); nit != nblist.end(); ++nit) {
        double J2 = log10((*nit)->getJeff2(state));

        MIN = (J2 < MIN) ? J2 : MIN;
        MAX = (J2 > MAX) ? J2 : MAX;

        J2s.push_back(J2);
    }

    // Prepare bins
    int BIN = ( (MAX-MIN)/_resolution_logJ2 + 0.5 ) + 1;
    vector< vector<double> > histJ2;
    histJ2.resize(BIN);

    // Execute binning
    vector< double > ::iterator jit;
    for (jit = J2s.begin(); jit < J2s.end(); ++jit) {

        int bin = int( (*jit-MIN)/_resolution_logJ2 + 0.5 );
        histJ2[bin].push_back(*jit);
    }

    vector< int > histN;
    histN.resize(BIN);
    for (int bin = 0; bin < BIN; ++bin) {
        histN[bin] = histJ2[bin].size();
    }
    FILE *out;
    string tag = boost::lexical_cast<string>(top->getDatabaseId())
               + "_INTEGRALS_" + ( (state == -1) ? "e" : "h" ) + ".dat";
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# IANALYZE: PAIR-INTEGRAL J2 HISTOGRAM\n");
    fprintf(out, "# STATE %1d\n", state);

    for (int bin = 0; bin < BIN; ++bin) {
        double J2 = MIN + bin*_resolution_logJ2;
        fprintf(out, "%4.7f %4d \n", J2, histN[bin]);
    }
    fclose(out);
}



}}



#endif
