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
    double _resolution_space;

};



void EAnalyze::Initialize(Topology *top, Property *opt) {

    string key = "options.eanalyze";

    _resolution_pairs = opt->get(key+".resolution_pairs").as< double >();
    _resolution_sites = opt->get(key+".resolution_sites").as< double >();
    _resolution_space = opt->get(key+".resolution_space").as< double >();
}

bool EAnalyze::EvaluateFrame(Topology *top) {

    // Calculate
    // ... Site-energy histogram, mean, width
    // ... Pair-energy histogram, mean, width
    // ... Site-energy correlation

    QMNBList &nblist = top->NBList();

    int state = -1;

    if (!top->Segments().size()) {
        cout << endl << "... ... Charge state " << state;
        cout << endl << "... ... ... No segments in topology. Skip ... "
             << flush;
    }
    else {
        SiteHist(top, state);
        SiteCorr(top, state);
    }

    

    if (!nblist.size()) {
        cout << endl << "... ... Charge state " << state;
        cout << endl << "... ... ... No pairs in topology. Skip ... "
             << flush;
    }
    else {
        PairHist(top, state);
    }

    

    state = +1;

    if (!top->Segments().size()) {
        cout << endl << "... ... Charge state " << state;
        cout << endl << "... ... ... No segments in topology. Skip ... "
             << flush;
    }
    else {
        SiteHist(top, state);
        SiteCorr(top, state);
    }

    if (!nblist.size()) {
        cout << endl << "... ... Charge state " << state;
        cout << endl << "... ... ... No pairs in topology. Skip ... "
             << flush;
    }
    else {
        PairHist(top, state);
    }

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
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# EANALYZE: SITE-ENERGY HISTOGRAM \n");
    fprintf(out, "# AVG %4.7f STD %4.7f MIN %4.7f MAX %4.7f \n", 
                    AVG,      STD,      MIN,      MAX);

    for (int bin = 0; bin < BIN; ++bin) {
        double E = MIN + bin*_resolution_sites;
        fprintf(out, "%4.7f %4d \n", E, histN[bin]);
    }
    fclose(out);
}


void EAnalyze::PairHist(Topology *top, int state) {

    QMNBList &nblist = top->NBList();
    QMNBList::iterator pit;

    double MIN = nblist.front()->Seg1()->getEMpoles(state)
               - nblist.front()->Seg2()->getEMpoles(state);
    double MAX = nblist.front()->Seg1()->getEMpoles(state)
               - nblist.front()->Seg2()->getEMpoles(state);
    double AVG = 0.0;
    double VAR = 0.0;
    double STD = 0.0;

    FILE *out_dEs;
    string tag_dEs = boost::lexical_cast<string>(top->getDatabaseId())
               + "_PAIRSLIST_" + ( (state == -1) ? "A" : "C" ) + ".dat";
    out_dEs = fopen(tag_dEs.c_str(), "w");

    // Collect site-energy differences from neighbourlist
    vector< double > dEs;    
    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {

        Segment *seg1 = (*pit)->Seg1();
        Segment *seg2 = (*pit)->Seg2();

        double dE = seg2->getEMpoles(state) - seg1->getEMpoles(state);

        MIN = (dE < MIN) ? dE : MIN;
        MIN = (-dE < MIN) ? -dE : MIN;
        MAX = (dE > MAX) ? dE : MAX;
        MAX = (-dE > MAX) ? -dE : MAX;
        AVG += dE / nblist.size();
        
        dEs.push_back(dE);
        dEs.push_back(-dE);

        fprintf(out_dEs, "%5d %5d %4.7f \n", seg1->getId(), seg2->getId(), dE);
    }
    fclose(out_dEs);
    
    // Prepare bins
    int BIN = int( (MAX-MIN)/_resolution_pairs + 0.5 ) + 1;
    vector< vector<double> > histE;
    histE.resize(BIN);
    
    // Execute binning, calc VAR
    vector< double > ::iterator eit;
    for (eit = dEs.begin(); eit < dEs.end(); ++eit) {

        int bin = int( ((*eit)-MIN)/_resolution_pairs + 0.5 );

        histE[bin].push_back((*eit));
        VAR += ((*eit) - AVG)*((*eit) - AVG) / nblist.size();
    }

    vector< int > histN;
    histN.resize(BIN);
    for (int bin = 0; bin < BIN; ++bin) {
        histN[bin] = histE[bin].size();
    }

    STD = sqrt(VAR);

    FILE *out;
    string tag = boost::lexical_cast<string>(top->getDatabaseId())
               + "_PAIRS_" + ( (state == -1) ? "A" : "C" ) + ".dat";
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# EANALYZE: PAIR-ENERGY HISTOGRAM \n");
    fprintf(out, "# AVG %4.7f STD %4.7f MIN %4.7f MAX %4.7f \n",
                    AVG,      STD,      MIN,      MAX);

    for (int bin = 0; bin < BIN; ++bin) {
        double E = MIN + bin*_resolution_pairs;
        fprintf(out, "%4.7f %4d \n", E, histN[bin]);
    }
    fclose(out);
}


void EAnalyze::SiteCorr(Topology *top, int state) {

    double AVG = 0.0;
    double VAR = 0.0;
    double STD = 0.0;

    vector< Segment* > ::iterator sit1;
    vector< Segment* > ::iterator sit2;    

    // Calculate mean site energy
    vector< double > Es;

    vector< Segment* > ::iterator sit;
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

        double E = (*sit)->getEMpoles(state);
        AVG += E / top->Segments().size();

        Es.push_back(E);
    }

    // Calculate variance
    vector< double > ::iterator eit;
    for (eit = Es.begin(); eit < Es.end(); ++eit) {

        VAR += ((*eit) - AVG)*((*eit) - AVG) / top->Segments().size();
    }

    // Collect inter-site distances, correlation product
    vector< double > Rs;
    vector< double > Cs;

    double MIN = +1e15;
    double MAX = -1e15;
    
    for (sit1 = top->Segments().begin(); sit1 < top->Segments().end(); ++sit1) {
    for (sit2 = sit1 + 1;                sit2 < top->Segments().end(); ++sit2) {

        double R = abs(top->PbShortestConnect((*sit1)->getPos(),
                                              (*sit2)->getPos()));

        MIN = (R < MIN) ? R : MIN;
        MAX = (R > MAX) ? R : MAX;

        double C = ((*sit1)->getEMpoles(state) - AVG)
                 * ((*sit2)->getEMpoles(state) - AVG);

        Rs.push_back(R);
        Cs.push_back(C);

    }}

    // Prepare bins
    int BIN = int( (MAX-MIN)/_resolution_space + 0.5 ) + 1;
    vector< vector<double> > histCs;
    histCs.resize(BIN);

    for (int i = 0; i < Rs.size(); ++i) {

        int bin = int((Rs[i] - MIN)/_resolution_space + 0.5);
        histCs[bin].push_back(Cs[i]);
    }

    // Calculate spatial correlation
    vector< double > histC;
    histC.resize(BIN);
    for (int bin = 0; bin < BIN; ++bin) {

        double corr = 0.0;
        for (int i = 0; i < histCs[bin].size(); ++i) {
            corr += histCs[bin][i] / VAR / histCs[bin].size();
        }
        histC[bin] = corr;
    }

    FILE *out;
    string tag = boost::lexical_cast<string>(top->getDatabaseId())
               + "_CORR_" + ( (state == -1) ? "A" : "C" ) + ".dat";
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# EANALYZE: SPATIAL SITE-ENERGY CORRELATION \n");
    fprintf(out, "# AVG %4.7f VAR %4.7f MIN_R %4.7f MAX_R %4.7f \n",
                    AVG,      VAR,      MIN,      MAX);

    for (int bin = 0; bin < BIN; ++bin) {
        double R = MIN + bin*_resolution_space;
        fprintf(out, "%4.7f %4.7f \n", R, histC[bin]);
    }
    fclose(out);
}







}}

#endif