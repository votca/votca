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
    void ParseEnergiesXML(Topology *top, Property *opt);
    void EvaluatePair(Topology *top, QMPair2 *pair);
    void CalculateRate(Topology *top, QMPair2 *pair, int state);


private:

    map<string, double> _seg_U_cC_nN_e;
    map<string, double> _seg_U_nC_nN_e;
    map<string, double> _seg_U_cN_cC_e;

    map<string, double> _seg_U_cC_nN_h;
    map<string, double> _seg_U_nC_nN_h;
    map<string, double> _seg_U_cN_cC_h;

    map<string, bool>   _seg_has_e;
    map<string, bool>   _seg_has_h;


    string _rateType;
    vec    _F;
    double _kT;
    double _omegaVib;
    int    _nMaxVib;

    int Factorial(int i);

};




void Rates2::Initialize(Topology *top, Property *options) {

    string key = "options.rates";

    /* ---- OPTIONS.XML Structure -----
     *
     * <rates>
     *
     *      <energiesXML>ENERGIES.XML</energiesXML>
     *
     *      <temperature></temperature>
     *      <field></field>
     *
     *      <method></method>
     *      <nmaxvib></nmaxvib>
     *      <omegavib></omegavib>
     *
     * </rates>
     *
     */

    // Control parameters
    double T = options->get(key+".temperature").as<double> ();
    _kT = 8.6173324e-5 * T;

    vec F = options->get(key+".field").as<vec>();
    _F = F;


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


    this->ParseEnergiesXML(top, options);

}


void Rates2::ParseEnergiesXML(Topology *top, Property *opt) {

    string key = "options.rates";
    string energiesXML = opt->get(key+".energiesXML").as<string> ();

    cout << endl
         << "... ... Site, reorg. energies from " << energiesXML << ". "
         << flush;

    Property alloc;
    load_property_from_xml(alloc, energiesXML.c_str());

    /* --- ENERGIES.XML Structure ---
     *
     * <topology>
     *
     *     <molecules>
     *          <molecule>
     *          <name></name>
     *
     *          <segments>
     *
     *              <segment>
     *              <name></name>
     *
     *              <!-- U_sG_sG, s->state, G->geometry !-->
     *
     *              <U_cC_nN_e></U_cC_nN_e>
     *              <U_cC_nN_h></U_cC_nN_h>
     *
     *              <U_nC_nN_e></U_nC_nN_e>
     *              <U_nC_nN_h></U_nC_nN_h>
     *
     *              <U_cN_cC_e></U_cN_cC_e>
     *              <U_cN_cC_h></U_cN_cC_h>
     *
     *              </segment>
     *
     *              <segment>
     *                  ...
     *
     */

    key = "topology.molecules.molecule";
    list<Property*> mols = alloc.Select(key);
    list<Property*> ::iterator molit;
    for (molit = mols.begin(); molit != mols.end(); ++molit) {

        key = "segments.segment";
        list<Property*> segs = (*molit)->Select(key);
        list<Property*> ::iterator segit;

        for (segit = segs.begin(); segit != segs.end(); ++segit) {

            string segName = (*segit)->get("name").as<string> ();

            bool has_e = false;
            bool has_h = false;

            double U_cC_nN_e = 0.0;
            double U_cC_nN_h = 0.0;
            double U_nC_nN_e = 0.0;
            double U_nC_nN_h = 0.0;
            double U_cN_cC_e = 0.0;
            double U_cN_cC_h = 0.0;

            if ( (*segit)->exists("U_cC_nN_e") &&
                 (*segit)->exists("U_nC_nN_e") &&
                 (*segit)->exists("U_cN_cC_e")    ) {

                U_cC_nN_e = (*segit)->get("U_cC_nN_e").as< double > ();
                U_nC_nN_e = (*segit)->get("U_nC_nN_e").as< double > ();
                U_cN_cC_e = (*segit)->get("U_cN_cC_e").as< double > ();

                has_e = true;
            }
            
            if ( (*segit)->exists("U_cC_nN_h") &&
                 (*segit)->exists("U_nC_nN_h") &&
                 (*segit)->exists("U_cN_cC_h")    ) {

                U_cC_nN_h = (*segit)->get("U_cC_nN_h").as< double > ();
                U_nC_nN_h = (*segit)->get("U_nC_nN_h").as< double > ();
                U_cN_cC_h = (*segit)->get("U_cN_cC_h").as< double > ();

                has_h = true;
            }

            _seg_U_cC_nN_e[segName] = U_cC_nN_e;
            _seg_U_nC_nN_e[segName] = U_nC_nN_e;
            _seg_U_cN_cC_e[segName] = U_cN_cC_e;
            _seg_has_e[segName] = has_e;

            _seg_U_cC_nN_h[segName] = U_cC_nN_h;
            _seg_U_nC_nN_h[segName] = U_nC_nN_h;
            _seg_U_cN_cC_h[segName] = U_cN_cC_h;
            _seg_has_h[segName] = has_h;
        }
    }
}


void Rates2::EvaluatePair(Topology *top, QMPair2 *qmpair) {

    cout << "\r... ... Evaluating pair " << qmpair->getId()+1 << flush;

    bool pair_has_e = false;
    bool pair_has_h = false;

    string segName1 = qmpair->first->getName();
    string segName2 = qmpair->second->getName();

    try {
        pair_has_e = _seg_has_e.at(segName1) && _seg_has_e.at(segName2);
        pair_has_h = _seg_has_h.at(segName1) && _seg_has_h.at(segName2);
    }
    catch (out_of_range) {
        cout << endl << "... ... WARNING: No energy information for pair ["
                     << segName1 << ", " << segName2 << "]. "
                     << "Skipping... " << endl;

        return;
    }

    if (pair_has_e) {
        this->CalculateRate(top, qmpair, -1);
    }
    if (pair_has_h) {
        this->CalculateRate(top, qmpair, +1);
    }
}


void Rates2::CalculateRate(Topology *top, QMPair2 *qmpair, int state) {

    const double NM2M    = 1.e-9;
    const double hbar_eV = 6.58211899e-16;

    Segment *seg1 = qmpair->first;
    Segment *seg2 = qmpair->second;

    double rate12 = 0.;                                       // 1->2
    double rate21 = 0.;                                       // 2->1

    double reorg12  = 0.0;                                    // 1->2
    double reorg21  = 0.0;                                    // 2->1
    double lOut     = qmpair->getLambdaO(state);              // 1->2 == + 2->1
    double dG_Site  = 0.0;                                    // 1->2 == - 2->1
    double dG_Field = - state * _F * qmpair->R() * NM2M;      // 1->2 == - 2->1
    double J2 = qmpair->getJeff2(state);                      // 1->2 == + 2->1

    if (state == -1) {

        // Charge hops Seg1 -> Seg2
        reorg12 = _seg_U_nC_nN_e[seg1->getName()]
                + _seg_U_cN_cC_e[seg2->getName()];

        // Charge hops Seg2 -> Seg1
        reorg21 = _seg_U_nC_nN_e[seg2->getName()]
                + _seg_U_cN_cC_e[seg1->getName()];

        // dG Seg1 -> Seg2
        dG_Site = _seg_U_cC_nN_e[seg2->getName()] + seg2->getEMpoles(state)
                - _seg_U_cC_nN_e[seg1->getName()] - seg1->getEMpoles(state);
    }

    else if (state == +1) {

        // Charge hops Seg1 -> Seg2
        reorg12 = _seg_U_nC_nN_h[seg1->getName()]
                + _seg_U_cN_cC_h[seg2->getName()];

        // Charge hops Seg2 -> Seg1
        reorg21 = _seg_U_nC_nN_h[seg2->getName()]
                + _seg_U_cN_cC_h[seg1->getName()];

        // dG Seg1 -> Seg2
        dG_Site = _seg_U_cC_nN_h[seg2->getName()] + seg2->getEMpoles(state)
                - _seg_U_cC_nN_h[seg1->getName()] - seg1->getEMpoles(state);
    }

    double dG = dG_Field + dG_Site;

    // +++++++++++++ //
    // JORTNER RATES //
    // +++++++++++++ //

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

    // ++++++++++++ //
    // MARCUS RATES //
    // ++++++++++++ //

    else if (_rateType == "marcus") {

        reorg12 = reorg12 + lOut;
        reorg21 = reorg21 + lOut;

        rate12 = J2 / hbar_eV * sqrt( M_PI / (reorg12*_kT) )
                * exp( - (+dG + reorg12)*(+dG + reorg12) / (4*_kT*reorg12) );

        rate21 = J2 / hbar_eV * sqrt( M_PI / (reorg21*_kT) )
                * exp( - (-dG + reorg21)*(-dG + reorg21) / (4*_kT*reorg21) );
    }

    qmpair->setRate12(rate12, state);
    qmpair->setRate21(rate21, state);
    qmpair->setIsPathCarrier(true, state);

}

int Rates2::Factorial(int i) {
    int k,j;
    k = 1;
    for (j=1; j<i; j++) { k = k*(j+1); }
    return k;
}

}}










#endif
