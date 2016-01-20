#ifndef _VELOCITY_H
#define _VELOCITY_H


#include <cstdlib>
#include <votca/xtp/qmcalculator.h>


namespace votca { namespace xtp {


class Velocity : public QMCalculator
{
public:

    string Identify() { return "Velocity"; }
    void Initialize(Property *opt);
    bool EvaluateFrame(Topology *top);

private:

    bool                _useCanonical;
    string              _includeSegsFile;
    vector<int>         _states;
    vector<Segment*>    _sites;

    double              _kT;

};



void Velocity::Initialize(Property *opt) {
    
    string key = "options.velocity";

    _states             = opt->get(key+".states").as< vector<int> >();
    _includeSegsFile    = opt->get(key+".include_sites").as< string >();
    _useCanonical       = (opt->get(key+".use_boltz_occ").as< int >() == 1)
                        ? true : false ;
    double T            = opt->get(key+".temperature").as< double >();

    if (_useCanonical) {
        cout << endl << "... ... Using canonical site occupations. " << flush;
    }

    cout << endl << "... ... Restricting sites as in file " << _includeSegsFile;
    
    _kT = 8.6173324e-5 * T;
}

bool Velocity::EvaluateFrame(Topology *top) {

    // ============================================== //
    // READ IN ALL SITES WHICH TO INCLUDE IN ANALYSIS //

    std::string line;
    std::ifstream intt;
    intt.open(_includeSegsFile.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {
            std::getline(intt, line);

            vector< string > split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);
            if ( !split.size()    ||
                  split[0] == "#" ||
                  split[0].substr(0,1) == "#" ) { continue; }

            int segId = boost::lexical_cast<int>(split[0]);
            _sites.push_back(top->getSegment(segId));
        }
    }
    else {
        throw std::runtime_error("No such file: '"+_includeSegsFile+"'.");
    }

    cout << endl << "... ... Including " << _sites.size() << " sites." << flush;
    

    for (unsigned int i = 0; i < _states.size(); ++i) {

        int state = _states[i];

        // ================================= //
        // CALCULATE SITE OCC. PROBABILITIES //

        double Z = 0;
        vector<Segment*> ::iterator sit;
        for (sit = _sites.begin(); sit < _sites.end(); ++sit) {

            Z += exp( - (*sit)->getSiteEnergy(state) / _kT );
        }

        double check_sum = 0.0;

        for (sit = _sites.begin(); sit < _sites.end(); ++sit) {

            double p = exp( - (*sit)->getSiteEnergy(state) / _kT ) / Z;
            check_sum += p;

            (*sit)->setOcc(p, state);
        }

        cout << endl << "1 = " << check_sum << flush;


        // ================================= //
        // CALCULATE AVG VELOCITY            //

        QMNBList &nblist = top->NBList();

        vector<Segment*> ::iterator sit1;
        vector<Segment*> ::iterator sit2;

        vec avg_v = vec(0,0,0);

        for (sit1 = _sites.begin(); sit1 < _sites.end(); ++sit1) {
        for (sit2 = sit1 + 1; sit2 < _sites.end(); ++sit2) {

          QMPair *qmp = nblist.FindPair((*sit1),(*sit2));

          if (qmp == NULL) {
              continue;
          }

          avg_v +=   (*sit1)->getOcc(state) * qmp->getRate12(state) * qmp->R();
          avg_v += - (*sit2)->getOcc(state) * qmp->getRate21(state) * qmp->R();
        }}


        avg_v *= 1e-9;

        cout << endl 
             << "... ... State " << state << ": v = " << avg_v
             << flush;

    }

    _sites.clear();
    
    return true;
}






}}



#endif
