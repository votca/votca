/*
 *            Copyright 2009-2016 The VOTCA Development Team
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

#ifndef VOTCA_XTP_IANALYZE_H
#define VOTCA_XTP_IANALYZE_H

#include <votca/xtp/qmcalculator.h>
#include <math.h>
#include <votca/xtp/qmpair.h>

namespace votca { namespace xtp {

class IAnalyze : public QMCalculator
{
public:

    string  Identify() { return "ianalyze"; }

    void    Initialize(Property *options);
    bool    EvaluateFrame(Topology *top);
    void    IHist(Topology *top, int state);
    void    IRdependence(Topology *top, int state);

private:

    double      _resolution_logJ2;
    vector<int> _states;
    double      _resolution_space;
    vector<QMPair::PairType> _pairtype;
    bool        _do_pairtype;
    bool        _do_IRdependence;

};


void IAnalyze::Initialize(Property *opt) {
    _do_pairtype=false;
    _do_IRdependence=false;
    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( opt );
    string key = "options." + Identify();
    _states = opt->get(key+".states").as< vector<int> >();
    _resolution_logJ2 = opt->get(key+".resolution_logJ2").as< double >();
    
    if ( opt->exists(key+".pairtype")) {
        _do_pairtype=true;
        string _store_string = opt->get(key+".pairtype").as<string> ();
        if (_store_string.find("Hopping") != std::string::npos) _pairtype.push_back(QMPair::Hopping);
        if (_store_string.find("SuperExchange") != std::string::npos) _pairtype.push_back(QMPair::SuperExchange);
        if (_store_string.find("SuperExchangeAndHopping") != std::string::npos) _pairtype.push_back(QMPair::SuperExchangeAndHopping);
        if (_store_string.find("Excitoncl") != std::string::npos) _pairtype.push_back(QMPair::Excitoncl);
        if (!_pairtype.size()){
            cout << endl << "... ... No pairtypes recognized will output all pairs. ";
             _do_pairtype=false;
        }
        //cout <<_pairtype.size()<<endl;
    }
    if ( opt->exists(key+".resolution_space")) {
        
        _resolution_space = opt->get(key+".resolution_space").as< double >();
        if (_resolution_space!=0.0) _do_IRdependence=true;
    }
   
}


bool IAnalyze::EvaluateFrame(Topology *top) {

    QMNBList &nblist = top->NBList();

    if (!nblist.size()) {
        cout << endl << "... ... No pairs in topology. Skip...";
        return 0;
    }
    
    if (_do_pairtype){
        bool pairs_exist=false;
        QMNBList::iterator nit;
        for (nit = nblist.begin(); nit != nblist.end(); ++nit) {
            QMPair::PairType pairtype=(*nit)->getType();
            if(std::find(_pairtype.begin(), _pairtype.end(), pairtype) != _pairtype.end()) {
                pairs_exist=true;
                break;
            }
        }
        if (!pairs_exist) {
        cout << endl << "... ... No pairs of given pairtypes in topology. Skip...";
        return 0;
    }
    }

    for (unsigned int i = 0; i < _states.size(); ++i) {
        this->IHist(top, _states[i]);
        if (_do_IRdependence){
        this->IRdependence(top, _states[i]);
        }
    }
    
    return true;
}


void IAnalyze::IHist(Topology *top, int state) {

    QMNBList &nblist = top->NBList();
    QMNBList::iterator nit;
   

    double MIN = std::numeric_limits<double>::max();
    double MAX = 0.0;
    
    // Collect J2s from pairs
    vector< double > J2s;
    //J2s.reserve(nblist.size()); //does not make a difference 
   
    for (nit = nblist.begin(); nit != nblist.end(); ++nit) {
        if(_do_pairtype){
            QMPair::PairType pairtype=(*nit)->getType();
            if(!(std::find(_pairtype.begin(), _pairtype.end(), pairtype) != _pairtype.end())){
                continue;
            }
        }
        double test = (*nit)->getJeff2(state);
        double J2;
 
        if(test <= 0) {continue;} // avoid -inf in output
        J2=log10(test);
    
        
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
    string name;
    if (state==-1) name="e";
    else if (state==1) name="h";
    else if (state==2) name="s";
    else if (state==3) name="t";
    string tag = boost::lexical_cast<string>("ianalyze.ihist_") +name + ".out";
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# IANALYZE: PAIR-INTEGRAL J2 HISTOGRAM\n");
    fprintf(out, "# STATE %1d\n", state);

    for (int bin = 0; bin < BIN; ++bin) {
        double J2 = MIN + bin*_resolution_logJ2;
        fprintf(out, "%4.7f %4d \n", J2, histN[bin]);
    }
    fclose(out);
}

void IAnalyze::IRdependence(Topology *top, int state) {
    
    QMNBList &nblist = top->NBList();
    QMNBList::iterator nit;

    double MIN  = log10(nblist.front()->getJeff2(state));
    double MAX  = log10(nblist.front()->getJeff2(state));
    double MINR = abs(nblist.front()->getR());
    double MAXR = abs(nblist.front()->getR());
    
    // Collect J2s from pairs
    vector< double > J2s;
    J2s.reserve(nblist.size());
    vector< double > distances;
    distances.reserve(nblist.size());

    for (nit = nblist.begin(); nit != nblist.end(); ++nit) {
        double J2 = log10((*nit)->getJeff2(state));

        MIN = (J2 < MIN) ? J2 : MIN;
        MAX = (J2 > MAX) ? J2 : MAX;
        
        double distance = abs((*nit)->getR());

        MINR = (distance < MINR) ? distance : MINR;
        MAXR = (distance > MAXR) ? distance : MAXR;
        
        distances.push_back(distance);
        J2s.push_back(J2);
    }
    
    // Prepare R bins
    int _pointsR = (MAXR-MINR)/_resolution_space;
    vector< vector<double> > rJ2;
    rJ2.resize(_pointsR);


    // Loop over distance
    for (int i = 0; i< _pointsR; ++i){
        double thisMINR = MINR + i*_resolution_space;
        double thisMAXR = MINR + (i+1)*_resolution_space;
        
        // now count Js that lie within this R range, calculate mean and sigma
        //double meanJ2 = 0;
        //double sigmaJ2 = 0;
        //int noJ2 = 0;
        
        vector< double > ::iterator jit;
        int j = 0;
        for (jit = J2s.begin(); jit < J2s.end(); ++jit) {
            if(thisMINR < distances[j] && distances[j] < thisMAXR){
                rJ2[i].push_back(*jit);  
            }
            j++;
        }
    }
    
    // make plot values
    vector< double > avgJ2;
    for (vector< vector<double> > ::iterator it = rJ2.begin() ; it != rJ2.end(); ++it){
        double thisavgJ2 = 0;
        for(unsigned i=0; i < (*it).size(); i++){
            thisavgJ2 += (*it)[i];
        }
        thisavgJ2 /= (*it).size();
        avgJ2.push_back(thisavgJ2);
    }
    vector< double > errJ2;
    int j = 0;
    for (vector< vector<double> > ::iterator it = rJ2.begin() ; it != rJ2.end(); ++it){
        double thiserrJ2 = 0;
        for(unsigned i=0; i < (*it).size(); i++){
            thiserrJ2 += ((*it)[i]-avgJ2[j])*((*it)[i]-avgJ2[j]);
        }
        thiserrJ2 /= (*it).size();
        thiserrJ2  = sqrt(thiserrJ2);
        errJ2.push_back(thiserrJ2);
        j++;
    }
    
    
    // print to file
    FILE *out;
    string name;
    if (state==-1) name="e";
    else if (state==1) name="h";
    else if (state==2) name="s";
    else if (state==3) name="t";
    string tag = boost::lexical_cast<string>("ianalyze.ispatial_") + name + ".out";
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# IANALYZE: SPATIAL DEPENDENCE OF log10(J2) [r,log10(J),error]\n");
    fprintf(out, "# STATE %1d\n", state);

    
    for (int i = 0; i < _pointsR; ++i) {
        double thisR = MINR + (i+0.5)*_resolution_space;
        fprintf(out, "%4.7f %4.7f %4.7f \n", thisR, avgJ2[i],  errJ2[i]);
    }
    fclose(out);
    

}



}}



#endif
