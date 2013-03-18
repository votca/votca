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

    vector<int> _states;

    double _site_avg;
    
    bool _skip_corr;
    bool _skip_sites;
    bool _skip_pairs;

};



void EAnalyze::Initialize(Topology *top, Property *opt) {

    string key = "options.eanalyze";

    _resolution_pairs = opt->get(key+".resolution_pairs").as< double >();
    _resolution_sites = opt->get(key+".resolution_sites").as< double >();
    _resolution_space = opt->get(key+".resolution_space").as< double >();

    if (opt->exists(key+".states")) {
        _states = opt->get(key+".states").as< vector<int> >();
    }
    else {
        _states.push_back(-1);
        _states.push_back(+1);
    }
    
    _skip_corr = opt->exists(key+".skip_correlation");
    _skip_sites = opt->exists(key+".skip_sites");
    _skip_pairs = opt->exists(key+".skip_pairs");

}

bool EAnalyze::EvaluateFrame(Topology *top) {

    // Calculate
    // ... Site-energy histogram, mean, width
    // ... Pair-energy histogram, mean, width
    // ... Site-energy correlation

    QMNBList &nblist = top->NBList();

    for (int i = 0; i < _states.size(); ++i) {

        int state = _states[i];
        cout << endl << "... ... Charge state " << state << flush;

        if (!top->Segments().size()) {            
            cout << endl << "... ... ... No segments in topology. Skip ... "
                 << flush;
        }
        else {
            // Site-energy histogram <> DOS
            if (_skip_sites) {
                cout << endl << "... ... ... Skip site-energy hist." << flush;
            }
            else {
                SiteHist(top, state);
            }
            
            // Site-energy correlation function
            if (_skip_corr) {
                cout << endl << "... ... ... Skip correlation ..." << flush;
            }
            else {
                SiteCorr(top, state);
            }
        }

        if (!nblist.size()) {
            cout << endl << "... ... ... No pairs in topology. Skip ... "
                 << flush;
        }
        else {
            // Site-energy-difference histogram <> Pair DOS
            if (_skip_pairs) {
                cout << endl << "... ... ... Skip pair-energy hist." << flush;
            }
            else {
                PairHist(top, state);
            }
        }
    }
}

void EAnalyze::SiteHist(Topology *top, int state) {

    vector< double > Es;
    Es.reserve(top->Segments().size());

    double MIN = top->Segments()[0]->getSiteEnergy(state);
    double MAX = top->Segments()[0]->getSiteEnergy(state);
    double AVG = 0.0;
    double VAR = 0.0;
    double STD = 0.0;

    // Collect energies from segments, calc AVG
    vector< Segment* > ::iterator sit;
    for (sit = top->Segments().begin(); 
         sit < top->Segments().end();
         ++sit) {

        double E = (*sit)->getSiteEnergy(state);

        MIN = (E < MIN) ? E : MIN;
        MAX = (E > MAX) ? E : MAX;
        AVG += E / top->Segments().size();
        
        Es.push_back(E);
    }

    _site_avg = AVG;

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
    
    
    tag = "x_y_z_e";
    out = fopen(tag.c_str(), "w");
    
    for (sit = top->Segments().begin(); 
         sit < top->Segments().end();
         ++sit) {
        
        if ((*sit)->getId() > 1792) { continue; }
        double E = (*sit)->getSiteEnergy(state);
        
        vector< Atom* > ::iterator ait;
        for (ait = (*sit)->Atoms().begin();
             ait < (*sit)->Atoms().end();
             ++ait) {
            
            Atom *atm = *ait;
            
            fprintf(out, "%3s %4.7f %4.7f %4.7f %4.7f\n",
                          (*sit)->getName().c_str(),
                          atm->getPos().getX(),
                          atm->getPos().getY(),
                          atm->getPos().getZ(),
                          E);            
        }
    }
    
    fclose(out);
    
}


void EAnalyze::PairHist(Topology *top, int state) {

    QMNBList &nblist = top->NBList();
    QMNBList::iterator pit;

    double MIN = nblist.front()->Seg1()->getSiteEnergy(state)
               - nblist.front()->Seg2()->getSiteEnergy(state);
    double MAX = nblist.front()->Seg1()->getSiteEnergy(state)
               - nblist.front()->Seg2()->getSiteEnergy(state);
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

        double dE = seg2->getSiteEnergy(state) - seg1->getSiteEnergy(state);

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

        double E = (*sit)->getSiteEnergy(state);
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

    vector< Fragment* > ::iterator fit1;
    vector< Fragment* > ::iterator fit2;

    cout << endl;

    for (sit1 = top->Segments().begin(); sit1 < top->Segments().end(); ++sit1) {

        cout << "\r... ... ..." << " Correlating segment ID = "
             << (*sit1)->getId() << flush;

    for (sit2 = sit1 + 1;                sit2 < top->Segments().end(); ++sit2) {

        double R = abs(top->PbShortestConnect((*sit1)->getPos(),
                                              (*sit2)->getPos()));

        for (fit1 = (*sit1)->Fragments().begin();
             fit1 < (*sit1)->Fragments().end();
             ++fit1) {
        for (fit2 = (*sit2)->Fragments().begin();
             fit2 < (*sit2)->Fragments().end();
             ++fit2) {

            double R_FF = abs(top->PbShortestConnect((*fit1)->getPos(),
                                                     (*fit2)->getPos()));

            if (R_FF < R) { R = R_FF; }

        }}


        MIN = (R < MIN) ? R : MIN;
        MAX = (R > MAX) ? R : MAX;

        double C = ((*sit1)->getSiteEnergy(state) - AVG)
                 * ((*sit2)->getSiteEnergy(state) - AVG);

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
    vector< double > histC_error;
    histC.resize(BIN);
    histC_error.resize(BIN);
    for (int bin = 0; bin < BIN; ++bin) {

        double corr = 0.0;
        double dcorr2 = 0.0;
        for (int i = 0; i < histCs[bin].size(); ++i) {
            corr += histCs[bin][i] / VAR;
            //corr2 += (histCs[bin][i] / VAR)*(histCs[bin][i] / VAR);
        }
        corr  = corr / histCs[bin].size();
        //corr2 = corr2 / histCs[bin].size();

        for (int i = 0; i < histCs[bin].size(); ++i) {
            dcorr2 += (histCs[bin][i]/VAR - corr)*(histCs[bin][i]/VAR - corr);
        }
        dcorr2 = dcorr2 / histCs[bin].size();
        
        histC[bin] = corr;
        histC_error[bin] = sqrt(dcorr2);
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
        fprintf(out, "%4.7f %4.7f %4.7f\n", R, histC[bin], histC_error[bin]);
    }
    fclose(out);
}







}}

#endif