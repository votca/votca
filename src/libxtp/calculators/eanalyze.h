/* 
 *            Copyright 2009-2017 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _VOTCA_XTP_EANALYZE_H
#define _VOTCA_XTP_EANALYZE_H

#include <votca/ctp/qmcalculator.h>
#include <math.h>
#include <votca/tools/tokenizer.h>


namespace votca { namespace xtp {

class EAnalyze : public ctp::QMCalculator
{
public:

    EAnalyze() { };
   ~EAnalyze() { };

    std::string Identify() { return "xeanalyze"; }

    void Initialize(tools::Property *opt);
    bool EvaluateFrame(ctp::Topology *top);
    void SiteHist(ctp::Topology *top, int state);
    void PairHist(ctp::Topology *top, int state);
    void SiteCorr(ctp::Topology *top, int state);

private:

    double _resolution_pairs;
    double _resolution_sites;
    double _resolution_space;
    std::string _distancemode;

    std::vector<int> _states;

    double _site_avg;
    
    bool _skip_corr;
    bool _skip_sites;
    bool _skip_pairs;
    
    bool _do_atomic_xyze;
    int  _atomic_first;
    int  _atomic_last;
    
    std::string _seg_pattern;
    std::vector<ctp::Segment*> _seg_shortlist;

};



void EAnalyze::Initialize( tools::Property *opt ) {
    _skip_corr=false;
    _skip_sites=false;
    _skip_pairs=false;
    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( opt, "xtp" );
    std::string key = "options." + Identify();
    if (opt->exists(key+".resolution_pairs")) {
    _resolution_pairs = opt->get(key+".resolution_pairs").as< double >();
    }
    else{
    _skip_pairs=true;
    }
    if (opt->exists(key+".resolution_sites")) {
    _resolution_sites = opt->get(key+".resolution_sites").as< double >();
    }
    else {
        _skip_sites=true;
    }
    if (opt->exists(key+".resolution_space")) {
    _resolution_space = opt->get(key+".resolution_space").as< double >();
    }
    else{
        _skip_corr=true;
    }

    if (opt->exists(key+".pattern")) {
        _seg_pattern = opt->get(key+".pattern").as<std::string>();
    }
    else _seg_pattern = "*";
    
    if (opt->exists(key+".states")) {
        _states = opt->get(key+".states").as< std::vector<int> >();
    }
    else {
        _states.push_back(-1);
        _states.push_back(+1);
    }
    
    if (opt->exists(key+".do_atomic_xyze")) {
        int do_xyze = opt->get(key+".do_atomic_xyze").as< int >();
        _do_atomic_xyze = (do_xyze == 1) ? true : false;
        _atomic_first = opt->get(key+".atomic_first").as< int >();
        _atomic_last  = opt->get(key+".atomic_last").as< int >();
    }
    else {
        _do_atomic_xyze = false;
    }
    
    if (opt->exists(key+".distancemode")) {
        // distancemode = segment / centreofmass
        _distancemode = opt->get(key+".distancemode").as< std::string >();
    }
    else{
         _distancemode = "segment";
    }
    if(_distancemode != "segment" && _distancemode != "centreofmass"){
        std::cout << "WARNING: distancemode has to be set to either 'segment' or to 'centreofmass'. Setting it to 'segment' now." << std::endl;
        _distancemode = "segment";
    }
    
    //_skip_corr = opt->exists(key+".skip_correlation");
    //_skip_sites = opt->exists(key+".skip_sites");
    //_skip_pairs = opt->exists(key+".skip_pairs");

}

bool EAnalyze::EvaluateFrame(ctp::Topology *top) {
    
    // Short-list segments according to pattern
    std::vector<ctp::Segment*>::iterator sit;
    for (sit=top->Segments().begin(); sit!=top->Segments().end(); ++sit) {
        std::string seg_name = (*sit)->getName();
        if (votca::tools::wildcmp(_seg_pattern.c_str(), seg_name.c_str())) {
            _seg_shortlist.push_back(*sit);
        }
    }
    std::cout << std::endl << "... ... Short-listed " << _seg_shortlist.size() 
         << " segments (pattern='" << _seg_pattern << "')" << std::flush;
    std::cout << std::endl << "... ... ... NOTE Statistics of site energies and spatial"
         << " correlations thereof are based on the short-listed segments only. "
         << std::flush;
    std::cout << std::endl << "... ... ...      "
         << "Statistics of site-energy differences operate on the full list." 
         << std::flush;

    // Calculate
    // ... Site-energy histogram, mean, width
    // ... Pair-energy histogram, mean, width
    // ... Site-energy correlation

    ctp::QMNBList &nblist = top->NBList();

    for (unsigned i = 0; i < _states.size(); ++i) {

        int state = _states[i];
        std::cout << std::endl << "... ... Charge state " << state << std::flush;

        if (!_seg_shortlist.size()) {            
            std::cout << std::endl << "... ... ... No segments short-listed. Skip ... "
                 << std::flush;
        }
        else {
            // Site-energy histogram <> DOS
            if (_skip_sites) {
                std::cout << std::endl << "... ... ... Skip site-energy hist." << std::flush;
            }
            else {
                SiteHist(top, state);
            }
            
            // Site-energy correlation function
            if (_skip_corr) {
                std::cout << std::endl << "... ... ... Skip correlation ..." << std::flush;
            }
            else {
                SiteCorr(top, state);
            }
        }

        if (!nblist.size()) {
            std::cout << std::endl << "... ... ... No pairs in topology. Skip ... "
                 << std::flush;
        }
        else {
            // Site-energy-difference histogram <> Pair DOS
            if (_skip_pairs) {
                std::cout << std::endl << "... ... ... Skip pair-energy hist." << std::flush;
            }
            else {
                PairHist(top, state);
            }
        }
    }
    
    return true;
}

void EAnalyze::SiteHist(ctp::Topology *top, int state) {

    std::vector< double > Es;
    Es.reserve(_seg_shortlist.size());

    double MIN = _seg_shortlist[0]->getSiteEnergy(state);
    double MAX = _seg_shortlist[0]->getSiteEnergy(state);
    double AVG = 0.0;
    double VAR = 0.0;
    double STD = 0.0;

    // Collect energies from segments, calc AVG
    std::vector< ctp::Segment* > ::iterator sit;
    for (sit = _seg_shortlist.begin(); 
         sit < _seg_shortlist.end();
         ++sit) {

        double E = (*sit)->getSiteEnergy(state);

        MIN = (E < MIN) ? E : MIN;
        MAX = (E > MAX) ? E : MAX;
        AVG += E / _seg_shortlist.size();
        
        Es.push_back(E);
    }

    _site_avg = AVG;

    // Prepare bins
    int BIN = int( (MAX-MIN)/_resolution_sites + 0.5 ) + 1;
    std::vector< std::vector<double> > histE;
    histE.resize(BIN);

    // Execute binning, calc VAR
    std::vector< double > ::iterator eit;
    for (eit = Es.begin(); eit < Es.end(); ++eit) {

        int bin = int( (*eit-MIN)/_resolution_sites + 0.5 );
        histE[bin].push_back(*eit);
        VAR += ((*eit) - AVG)*((*eit) - AVG) / _seg_shortlist.size();
    }

    std::vector< int > histN;
    histN.resize(BIN);
    for (int bin = 0; bin < BIN; ++bin) {
        histN[bin] = histE[bin].size();
    }

    STD = sqrt(VAR);

    FILE *out;
    std::string statename;
    if (state==-1){
        statename="e";
    }
    else if (state==1){
        statename="h";
    }
    else if (state==2){
        statename="s";
    }
    else if (state==3){
        statename="t";
    }
    
    
    std::string tag = boost::lexical_cast<std::string>("eanalyze.sitehist_") +statename+ ".out";
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# EANALYZE: SITE-ENERGY HISTOGRAM \n");
    fprintf(out, "# AVG %4.7f STD %4.7f MIN %4.7f MAX %4.7f \n", 
                    AVG,      STD,      MIN,      MAX);

    for (int bin = 0; bin < BIN; ++bin) {
        double E = MIN + bin*_resolution_sites;
        fprintf(out, "%4.7f %4d \n", E, histN[bin]);
    }
    fclose(out);
    
    // Write "seg x y z energy" with atomic {x,y,z}
    if (_do_atomic_xyze) {
        tag = (state == -1) ? "eanalyze.landscape_e.out" : "eanalyze.landscape_h.out";
        out = fopen(tag.c_str(), "w");

        for (sit = _seg_shortlist.begin(); 
             sit < _seg_shortlist.end();
             ++sit) {

            if ((*sit)->getId() < _atomic_first) { continue; }
            if ((*sit)->getId() == _atomic_last) { break; }
            double E = (*sit)->getSiteEnergy(state);

            std::vector< ctp::Atom* > ::iterator ait;
            for (ait = (*sit)->Atoms().begin();
                 ait < (*sit)->Atoms().end();
                 ++ait) {

                ctp::Atom *atm = *ait;

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
}


void EAnalyze::PairHist(ctp::Topology *top, int state) {

    ctp::QMNBList &nblist = top->NBList();
    ctp::QMNBList::iterator pit;

    double MIN = nblist.front()->Seg1()->getSiteEnergy(state)
               - nblist.front()->Seg2()->getSiteEnergy(state);
    double MAX = nblist.front()->Seg1()->getSiteEnergy(state)
               - nblist.front()->Seg2()->getSiteEnergy(state);
    double AVG = 0.0;
    double VAR = 0.0;
    double STD = 0.0;

    FILE *out_dEs;
    std::string tag_dEs = boost::lexical_cast<std::string>("eanalyze.pairlist_") + ( (state == -1) ? "e" : "h" ) + ".out";
    out_dEs = fopen(tag_dEs.c_str(), "w");

    // Collect site-energy differences from neighbourlist
    std::vector< double > dEs;    
    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {

        ctp::Segment *seg1 = (*pit)->Seg1();
        ctp::Segment *seg2 = (*pit)->Seg2();

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
    std::vector< std::vector<double> > histE;
    histE.resize(BIN);
    
    // Execute binning, calc VAR
    std::vector< double > ::iterator eit;
    for (eit = dEs.begin(); eit < dEs.end(); ++eit) {

        int bin = int( ((*eit)-MIN)/_resolution_pairs + 0.5 );

        histE[bin].push_back((*eit));
        VAR += ((*eit) - AVG)*((*eit) - AVG) / nblist.size();
    }

    std::vector< int > histN;
    histN.resize(BIN);
    for (int bin = 0; bin < BIN; ++bin) {
        histN[bin] = histE[bin].size();
    }

    STD = sqrt(VAR);

    FILE *out;
    std::string statename;
    if (state==-1){
        statename="e";
    }
    else if (state==1){
        statename="h";
    }
    else if (state==2){
        statename="s";
    }
    else if (state==3){
        statename="t";
    }
    std::string tag = boost::lexical_cast<std::string>("eanalyze.pairhist_") + statename + ".out";
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


void EAnalyze::SiteCorr(ctp::Topology *top, int state) {

    double AVG = 0.0;
    double AVGESTATIC = 0.0;
    double VAR = 0.0;
    double STD = 0.0;

    std::vector< ctp::Segment* > ::iterator sit1;
    std::vector< ctp::Segment* > ::iterator sit2;    

    // Calculate mean site energy
    std::vector< double > Es;

    std::vector< ctp::Segment* > ::iterator sit;
    for (sit = _seg_shortlist.begin();
         sit < _seg_shortlist.end();
         ++sit) {

        double E = (*sit)->getSiteEnergy(state);
        AVG += E / _seg_shortlist.size();
        
        AVGESTATIC += (*sit)->getEMpoles(state) / top->Segments().size();

        Es.push_back(E);
    }
    
    // Calculate variance
    std::vector< double > ::iterator eit;
    for (eit = Es.begin(); eit < Es.end(); ++eit) {

        VAR += ((*eit) - AVG)*((*eit) - AVG) / _seg_shortlist.size();
    }
    
    STD = sqrt(VAR);

    // Collect inter-site distances, correlation product
    std::vector< double > Rs;
    std::vector< double > Cs;

    double MIN = +1e15;
    double MAX = -1e15;

    std::vector< ctp::Fragment* > ::iterator fit1;
    std::vector< ctp::Fragment* > ::iterator fit2;

    std::cout << std::endl;
    
    FILE *corr_out;
    std::string statename;
    if (state==-1){
        statename="e";
    }
    else if (state==1){
        statename="h";
    }
    else if (state==2){
        statename="s";
    }
    else if (state==3){
        statename="t";
    }
    std::string corrfile = boost::lexical_cast<std::string>("eanalyze.sitecorr.atomic_") + statename+ ".out";
    corr_out = fopen(corrfile.c_str(), "w");
    
    for (sit1 = _seg_shortlist.begin(); sit1 < _seg_shortlist.end(); ++sit1) {

        std::cout << "\r... ... ..." << " Correlating segment ID = "
             << (*sit1)->getId() << std::flush;

    for (sit2 = sit1 + 1;                sit2 < _seg_shortlist.end(); ++sit2) {

        double R = abs(top->PbShortestConnect((*sit1)->getPos(),
                                              (*sit2)->getPos()));

        if(_distancemode == "segment"){
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
        }
    

        MIN = (R < MIN) ? R : MIN;
        MAX = (R > MAX) ? R : MAX;

        double C = ((*sit1)->getSiteEnergy(state) - AVG)
                 * ((*sit2)->getSiteEnergy(state) - AVG);

        Rs.push_back(R);
        Cs.push_back(C);
        
        fprintf(corr_out, "%+1.7f %+1.7f\n", R, C);

    }}
    
    fclose(corr_out);

    // Prepare bins
    int BIN = int( (MAX-MIN)/_resolution_space + 0.5 ) + 1;
    std::vector< std::vector<double> > histCs;
    histCs.resize(BIN);

    for (unsigned i = 0; i < Rs.size(); ++i) {

        int bin = int((Rs[i] - MIN)/_resolution_space + 0.5);
        histCs[bin].push_back(Cs[i]);
    }

    // Calculate spatial correlation
    std::vector< double > histC;
    std::vector< double > histC_error;
    histC.resize(BIN);
    histC_error.resize(BIN);
    for (int bin = 0; bin < BIN; ++bin) {

        double corr = 0.0;
        double dcorr2 = 0.0;
        for (unsigned i = 0; i < histCs[bin].size(); ++i) {
            corr += histCs[bin][i] / VAR;
            //corr2 += (histCs[bin][i] / VAR)*(histCs[bin][i] / VAR);
        }

        corr  = corr / histCs[bin].size();
        //corr2 = corr2 / histCs[bin].size();

        for (unsigned i = 0; i < histCs[bin].size(); ++i) {
            dcorr2 += (histCs[bin][i]/VAR/histCs[bin].size() - corr)*(histCs[bin][i]/VAR/histCs[bin].size() - corr);
        }


        histC[bin] = corr;
        
        // error on mean value
        dcorr2 = dcorr2 / histCs[bin].size() / (histCs[bin].size()-1);
        histC_error[bin] = sqrt(dcorr2);
    }

    FILE *out;
    
    if (state==-1){
        statename="e";
    }
    else if (state==1){
        statename="h";
    }
    else if (state==2){
        statename="s";
    }
    else if (state==3){
        statename="t";
    }
    std::string tag = boost::lexical_cast<std::string>("eanalyze.sitecorr_") + statename + ".out";
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# EANALYZE: SPATIAL SITE-ENERGY CORRELATION \n");
    fprintf(out, "# AVG %4.7f STD %4.7f MIN_R %4.7f MAX_R %4.7f  AVGESTATIC %4.7f\n",
                    AVG,      STD,      MIN,      MAX,      AVGESTATIC);

    for (int bin = 0; bin < BIN; ++bin) {
        double R = MIN + bin*_resolution_space;
        fprintf(out, "%4.7f %4.7f %4.7f\n", R, histC[bin], histC_error[bin]);
    }
    fclose(out);
}

}}

#endif // _VOTCA_XTP_EANALYZE_H
