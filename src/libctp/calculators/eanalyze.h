#ifndef EANALYZE_H
#define EANALYZE_H

#include <votca/ctp/qmcalculator.h>
#include <math.h>


namespace votca { namespace ctp {

class EAnalyze : public QMCalculator
{
public:

    EAnalyze() { };
   ~EAnalyze() { };

    string Identify() { return "EAnalyze"; }

    void Initialize(Topology *top, Property *opt);
    bool EvaluateFrame(Topology *top);
    void SiteHist(Topology *top, int state);
    void PairHist(Topology *top, int state);
    void SiteCorr(Topology *top, int state);

private:

    double _resolution_pairs;
    double _resolution_sites;

};



void EAnalyze::Initialize(Topology *top, Property *opt) {

    string key = "options.eanalyze";

    _resolution_pairs = opt->get(key+".resolution_pairs").as< double >();
    _resolution_sites = opt->get(key+".resolution_sites").as< double >();

}

bool EAnalyze::EvaluateFrame(Topology *top) {

    // Calculate
    // ... Site-energy histogram, mean, width
    // ... Pair-energy histogram, mean, width
    // ... Site-energy correlation

    int state = -1;
    SiteHist(top, state);


    state = +1;
    SiteHist(top, state);

}

void EAnalyze::SiteHist(Topology *top, int state) {

    vector< double > Es;
    Es.reserve(top->Segments().size());

    double MIN = top->Segments()[0]->getEMpoles(state);
    double MAX = top->Segments()[0]->getEMpoles(state);
    double AVG = 0.0;
    double VAR = 0.0;
    double STD = 0.0;

    // Collect energies from segments, calc AVG
    vector< Segment* > ::iterator sit;
    for (sit = top->Segments().begin(); 
         sit < top->Segments().end();
         ++sit) {

        double E = (*sit)->getEMpoles(state);

        MIN = (E < MIN) ? E : MIN;
        MAX = (E > MAX) ? E : MAX;
        AVG += E / top->Segments().size();
        
        Es.push_back(E);
    }

    // Prepare bins
    int BIN = int( (MAX-MIN)/_resolution_sites + 0.5 ) + 1;
    vector< vector<double> > histE;
    histE.resize(BIN);

    // Execute binning, calc VAR
    vector< double > ::iterator eit;
    for (eit = Es.begin(); eit < Es.end(); ++eit) {

        int bin = int( (*eit-MIN)/_resolution_sites + 0.5 );
        histE[bin].push_back(*eit);
        VAR += ((*eit) - AVG)*((*eit) - AVG) / top->Segments().size();
    }

    vector< int > histN;
    histN.resize(BIN);
    for (int bin = 0; bin < BIN; ++bin) {
        histN[bin] = histE[bin].size();
    }

    STD = sqrt(VAR);

    FILE *out;
    string tag = boost::lexical_cast<string>(top->getDatabaseId())
               + "_SITES_" + ( (state == -1) ? "A" : "C" ) + ".dat";

    fprintf(out, "# EANALYZE: SITE-ENERGY HISTOGRAM \n");
    fprintf(out, "# AVG %4.7f STD %4.7f MIN %4.7f MAX %4.7f \n", 
                    AVG,      STD,      MIN,      MAX);

    for (int bin = 0; bin < BIN; ++bin) {
        double E = MIN + bin*_resolution_sites;
        fprintf(out, "%4.7f %4d \n", E, histN[bin]);
    }
}

void EAnalyze::PairHist(Topology *top, int state) {
    ;
}

void EAnalyze::SiteCorr(Topology *top, int state) {
    ;
}







}}

#endif