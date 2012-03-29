#ifndef RATES2_H
#define RATES2_H

#include <votca/ctp/paircalculator2.h>
#include <math.h>

namespace votca { namespace ctp {

class Rates2 : public PairCalculator2
{
public:

    Rates2() { };
   ~Rates2() { };

    string Identify() { return "Rates"; }

    void Initialize(Topology *top, Property *options);
    void EvaluatePair(Topology *top, QMPair2 *pair);


private:

    string _rateType;
    vec    _E;
    double _kT;
    double _omegaVib;
    int    _nMaxVib;

    int Factorial(int i);

};




void Rates2::Initialize(Topology *top, Property *options) {

    string key = "options.rates";


    // Control parameters
    double T = options->get(key+".temperature").as<double> ();
    _kT = 8.6173324e-5 * T;

    vec E = options->get(key+".field").as<vec>();
    _E = E;


    // Method
    if (options->exists(key+".method")) {
        _rateType = options->get(key+".method").as<string> ();
    }
    else {
        cout << endl 
             << "... ... ERROR: No method to calculate rates specified. "
             << endl;
        throw std::runtime_error("Missing input in options file.");
    }
    if (_rateType != "marcus" && _rateType != "jortner") {
        cout << endl
             << "... ... ERROR: Unknown rate type '" << _rateType << "' "
             << endl;
        throw std::runtime_error("Faulty input in options file.");
    }


    // Vibrational quanta
    if (options->exists(key+".nmaxvib")) {
        _nMaxVib = options->get(key+".nmaxvib").as<double> ();
    }
    else {
        _nMaxVib = 20;
        cout << endl << "... ... WARNING: No cut-off number for QM vibrations "
                        "provided, using default 20.";
    }
    if (options->exists(key+".omegavib")) {
        _omegaVib = options->get(key+".omegavib").as<double> ();
    }
    else {
        _omegaVib = 0.2;
        cout << endl << "... ... WARNING: No QM vibration frequency provided, "
                        "using default 0.2eV.";
    }
}

void Rates2::EvaluatePair(Topology *top, QMPair2 *qmpair) {

    cout << "\r... ... Evaluating pair " << qmpair->getId()+1 << flush;

    double NM2M = 1.e-9;
    const double hbar_eV = 6.58211899e-16;

    // TODO Decide how to include different carrier types
    int q = -1;


    Segment *seg1 = qmpair->first;
    Segment *seg2 = qmpair->second;

    double rate12 = 0.;
    double rate21 = 0.;

    double J2 = qmpair->calcJeff2();
    //cout << endl << "J2 " << J2;

    double reorg12 = seg1->getLambdaIntra(q, 0) + seg2->getLambdaIntra(0, q);
    double reorg21 = seg1->getLambdaIntra(0, q) + seg2->getLambdaIntra(q, 0);

    //cout << endl << "lambda " << reorg12 << " --- " << reorg21;

    
    double dG_Field = - _E * qmpair->R() * NM2M;
    //cout << endl << "dG Field" << dG_Field;
    double dG_Site = seg2->getESite(q) - seg1->getESite(q);
    //cout << endl << "dG Site " << dG_Site;
    double dG = dG_Field + dG_Site;

    double lOut = qmpair->getLambdaO();
    //cout << endl << "lambda Outer" << lOut;
    
    if (_rateType == "jortner" && lOut < 0.) {
        cout << endl
             << "... ... ERROR: Pair " << qmpair->getId() << " has negative "
                "outer-sphere reorganization energy. Cannot calculate Jortner "
                "rates. "
             << endl;
        throw std::runtime_error("");
    }
    else if (_rateType == "jortner" && lOut < 0.01) {
        cout << endl
             << "... ... WARNING: Pair " << qmpair->getId() << " has small "
                "outer-sphere reorganization energy (" << lOut << "eV). Could "
                "lead to over-estimated Jortner rates."
             << endl;
    }

    if (_rateType == "jortner") {

        double huang_rhys12 = reorg12 / _omegaVib;
        double huang_rhys21 = reorg21 / _omegaVib;

        int nvib;
        for (nvib = 0; nvib <= _nMaxVib; nvib++) {

            // Hopping from Seg1 -> Seg2
            rate12 += 1 / hbar_eV * sqrt( M_PI / (lOut*_kT) )
                    * J2 * exp(-huang_rhys12) * pow(huang_rhys12, nvib) /
                           Factorial(nvib)
                    * exp( -pow( (dG + nvib*_omegaVib + lOut) , 2 ) /
                           (4*_kT*lOut) );

            // Hopping from Seg2 -> Seg1
            rate21 += 1 / hbar_eV * sqrt( M_PI / (lOut*_kT) )
                    * J2 * exp(-huang_rhys21) * pow(huang_rhys21, nvib) /
                           Factorial(nvib)
                    * exp( -pow( (-dG + nvib*_omegaVib + lOut) , 2 ) /
                           (4*_kT*lOut) );            
        }
    }



    if (_rateType == "marcus") {

        reorg12 = reorg12 + lOut;
        reorg21 = reorg21 + lOut;

        rate12 = 1 / hbar_eV * sqrt( M_PI / (reorg12*_kT) )
                * J2 * exp( - pow( (dG + reorg12) , 2) /
                              (4*_kT*reorg12) );
        rate21 = 1 / hbar_eV * sqrt( M_PI / (reorg21*_kT) )
                * J2 * exp( - pow( (-dG + reorg21) , 2) /
                              (4*_kT*reorg21) );
    }

    //cout << endl << "Rate12 " << rate12;
    //cout << endl << "Rate21 " << rate21;

    qmpair->setRate12(q, rate12);
    qmpair->setRate21(q, rate21); 

}

int Rates2::Factorial(int i) {
    int k,j;
    k=1;
    for (j=1;j<i;j++) { k=k*(j+1); }
    return k;
}

}}










#endif
