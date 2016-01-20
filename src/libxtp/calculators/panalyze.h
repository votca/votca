#ifndef _VOTCA_XTP_PANALYZE_H
#define _VOTCA_XTP_PANALYZE_H

#include <votca/tools/globals.h>
#include <votca/xtp/qmcalculator.h>
#include <votca/xtp/qmpair.h>

#include <math.h>


namespace votca { namespace xtp {

class PAnalyze : public QMCalculator
{
public:

    PAnalyze() { };
   ~PAnalyze() { };

    string Identify() { return "panalyze"; }

    void Initialize(Property *opt);
    bool EvaluateFrame(Topology *top);
    void SiteConnection(Topology *top);
    void CoordinationNumber(Topology *top);

private:

    double _resolution_space;

};



void PAnalyze::Initialize( Property *opt ) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( opt );
    string key = "options." + Identify();

    _resolution_space = opt->get(key+".resolution_space").as< double >();
}

bool PAnalyze::EvaluateFrame(Topology *top) {

    // Calculate
    // ... Number of neighbors distribution, mean, width
    // ... Distance dependent probability for sites being connected
    

    QMNBList &nblist = top->NBList();
    if (!nblist.size()) {
        cout << endl << "... ... No pairs in topology. Skip...";
        return 0;
    }
    cout << endl << "... ... Evaluating coordination number (=number of neighbors) distribution." << endl;
    CoordinationNumber(top);
    cout << endl << "... ... Evaluating centre of mass distance-dependent probability of two sites being connected." << endl;
    SiteConnection(top);
    
    return true;
}


void PAnalyze::CoordinationNumber(Topology *top) {
    vector <Segment*> segments = top->Segments();
    int numberofmolecules = segments.size();
    vector<int> numberofneighbors;
    numberofneighbors.resize(numberofmolecules);
    cout << "... ... " << numberofmolecules << " sites to be analyzed."<< endl;
    
    QMNBList &nblist = top->NBList();
    QMNBList::iterator nit;
    for (nit = nblist.begin(); nit != nblist.end(); ++nit) {
        int seg1= (*nit)->Seg1()->getId();
        int seg2= (*nit)->Seg2()->getId();
        if(seg1-1 > numberofmolecules || seg2-1 > numberofmolecules){
            cout << "WARNING: state file does not contain consecutive segment enumeration up to " << numberofmolecules << ". Something is wrong." << endl;
            break;
        }
        numberofneighbors[seg1]++;
        numberofneighbors[seg2]++;
    }
    
    int MIN = numberofmolecules;
    int MAX = 0;
    double AVG = 0;
    for(int i=0; i<numberofmolecules; i++){
        MIN = (numberofneighbors[i] < MIN) ? numberofneighbors[i] : MIN;
        MAX = (numberofneighbors[i] > MAX) ? numberofneighbors[i] : MAX;
        AVG += numberofneighbors[i];
    }
    AVG /= numberofmolecules;
    cout << "... ... Coordination numbers between " << MIN << " and " << MAX << " found. Average: " << AVG << " neighbors." << endl;
    
    // make histogram
    vector<int> frequencies;
    for(int i=MIN; i<=MAX; i++) {
        int frequency = 0;
        for(int mol=0; mol<numberofmolecules; mol++){
            if(numberofneighbors[mol] == i){
                frequency++;
            }
        }
        frequencies.push_back(frequency);
    }
    
    // calculate standard deviation
    double STD = 0;
    for(int mol=0; mol<numberofmolecules; mol++){
        STD += (numberofneighbors[mol]-AVG)*(numberofneighbors[mol]-AVG);
    }
    STD /= numberofmolecules;
    STD  = sqrt(STD);
    
    
    FILE *out;
    string tag = boost::lexical_cast<string>("panalyze.coordination.out");
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# PANALYZE: NUMBER OF NEIGHBORS. \n");
    fprintf(out, "# AVG %4.7f STD %4.7f MIN %5d MAX %5d \n", 
                    AVG,      STD,      MIN,      MAX);
    for (unsigned i=0; i<frequencies.size(); i++) {
        fprintf(out, "%5d %5d \n", MIN+i, frequencies[i]);
    }
    fclose(out);
    
}


void PAnalyze::SiteConnection(Topology *top) {
    // get neighbours from topology
    QMNBList &nblist = top->NBList();
    QMNBList::iterator nit;

    vector< vector<int> > segids;
    vector<double> segdist;
    segids.reserve(nblist.size());
    segdist.reserve(nblist.size());
    
    for (nit = nblist.begin(); nit != nblist.end(); ++nit) {
        int seg1= (*nit)->Seg1()->getId();
        int seg2= (*nit)->Seg2()->getId();
        vector<int> thissegids;
        thissegids.push_back(seg1);
        thissegids.push_back(seg2);
        segids.push_back(thissegids);
        
        double thissegdist = abs((*nit)->R());
        segdist.push_back(thissegdist);
    }
    

    
    // get also non-neighbour connections from topology
    vector <Segment*> segments = top->Segments();
    vector <Segment*>::iterator seg1;
    vector <Segment*>::iterator seg2;
    
    // get MINR and MAXR
    double MINR = abs(nblist.front()->getR());
    double MAXR = abs(nblist.front()->getR());
    for (nit = nblist.begin(); nit != nblist.end(); ++nit) {
    double distance = abs((*nit)->getR());
    MINR = (distance < MINR) ? distance : MINR;
    MAXR = (distance > MAXR) ? distance : MAXR;
    }

    // Prepare R bins
    int _pointsR = ((MAXR-MINR)/_resolution_space+0.5) +1;
    cout << "... ... minimal R: " << MINR << " nm" << endl;
    cout << "... ... maximal R: " << MAXR << " nm" << endl;
    cout << "... ... R points:  " << _pointsR << endl;
    vector<double> Rconnected;
    vector<double> Rtotal;
    vector<double> Rprobability;
    Rconnected.resize(_pointsR);
    Rtotal.resize(_pointsR);
    Rprobability.resize(_pointsR);


    // find which sites are in the neighbor list and which are not    
    for (seg1 = segments.begin(); seg1!= segments.end(); seg1++){
        cout << "\r... ... ..." << " checking segment ID = "
             << (*seg1)->getId() << flush;

        for (seg2 = seg1; seg2!= segments.end(); seg2++){ // for (seg2 = segments.begin(); seg2!= segments.end(); seg2++):q
            vec r1 = (*seg1)->getPos();
            vec r2 = (*seg2)->getPos();
            double distance = abs( top->PbShortestConnect(r1, r2));
            if(MINR <= distance && distance <= MAXR){
                int inpairlist = 0;
                vector< vector<int> >::iterator thissegids;
                for(unsigned i = 0; i<segids.size(); i++){
                    if(((*seg1)->getId()==segids[i][0] && (*seg2)->getId() == segids[i][1]) || ((*seg1)->getId()==segids[i][1] && (*seg2)->getId() == segids[i][2])){
                        inpairlist = 1;
                        break;
                    }
                }
                int j = floor((distance-MINR)/_resolution_space);
                if(j < 0 || j > _pointsR) {cout << "WARNING: R bins seem incorrect. This should not happen. r = " << distance << ", j = " << j; break;}
                Rtotal[j] ++;
                if(inpairlist == 1){
                    Rconnected[j] ++;
                }
            }
        }
        
        
        
    }

    // evaluate connection probability function
    for (int j = 0; j< _pointsR; ++j){
        //double thisMINR = MINR + j*_resolution_space;
        if(Rconnected[j] == 0) {
            Rprobability[j]=0;
        }
        else{
            Rprobability[j] = Rconnected[j] / Rtotal[j];
        }
        // cout << thisMINR << "    " << Rprobability[j] << endl;
    }
    
    Rprobability.push_back(0);
    
    cout << endl << "... ... Done with evaluation. Now writing output files.";

    FILE *out;
    string tag = boost::lexical_cast<string>("panalyze.distanceprobability.out");
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# PANALYZE: CONNECTION PROBABILITY DEPENDING ON CENTRE-OF-MASS DISTANCE. \n");
    for (int j = 0; j < _pointsR+1; ++j) {
        double thisMINR = MINR + (j+0.5)*_resolution_space;
        fprintf(out, "%4.7f %4.7f \n", thisMINR, Rprobability[j]);
    }
    fclose(out);
}




}}

#endif // _VOTCA_XTP_PANALYZE_H