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


#ifndef __VOTCA_XTP_RATESCALC_H
#define __VOTCA_XTP_RATESCALC_H

#include <votca/xtp/paircalculator.h>
#include <cmath>
#include <complex>


namespace votca { namespace xtp {
   
class Rates : public PairCalculator
{
public:

    Rates() { };
   ~Rates() { };

    std::string Identify() { return "rates"; }

    void Initialize(tools::Property *options);
    void ParseEnergiesXML(Topology *top, tools::Property *opt);
    void EvaluatePair(Topology *top, QMPair *pair);
    void CalculateRate(Topology *top, QMPair *pair, int state);


private:

    std::map<std::string, double> _seg_U_cC_nN_e;
    std::map<std::string, double> _seg_U_nC_nN_e;
    std::map<std::string, double> _seg_U_cN_cC_e;

    std::map<std::string, double> _seg_U_cC_nN_h;
    std::map<std::string, double> _seg_U_nC_nN_h;
    std::map<std::string, double> _seg_U_cN_cC_h;

    std::map<std::string, bool>   _seg_has_e;
    std::map<std::string, bool>   _seg_has_h;


    std::string _rateType;
    tools::vec    _F;
    double _kT;
    double _omegaVib;
    int    _nMaxVib;
    double _kondo;

    int Factorial(int i);

};

std::complex<double> ccgamma(std::complex<double> z,int OPT)
{
    std::complex<double> g;
    double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
    double na=0,t,x1=0,y1,sr,si;
    //int i,j,k;
    int j,k;

    static double a[] = {
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590};

    x = real(z);
    y = imag(z);
    if (x > 171) return std::complex<double>(1e308,0);
    if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
        return std::complex<double>(1e308,0);
    else if (x < 0.0) {
        x1 = x;
        y1 = y;
        x = -x;
        y = -y;
    }
    x0 = x;
    if (x <= 7.0) {
        na = (int)(7.0-x);
        x0 = x+na;
    }
    q1 = sqrt(x0*x0+y*y);
    th = atan(y/x0);
    gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
    gi = th*(x0-0.5)+y*log(q1)-y;
    for (k=0;k<10;k++){
        t = pow(q1,-1.0-2.0*k);
        gr += (a[k]*t*cos((2.0*k+1.0)*th));
        gi -= (a[k]*t*sin((2.0*k+1.0)*th));
    }
    if (x <= 7.0) {
        gr1 = 0.0;
        gi1 = 0.0;
        for (j=0;j<na;j++) {
            gr1 += (0.5*log((x+j)*(x+j)+y*y));
            gi1 += atan(y/(x+j));
        }
        gr -= gr1;
        gi -= gi1;
    }
    if (x1 <= 0.0) {
        q1 = sqrt(x*x+y*y);
        th1 = atan(y/x);
        sr = -sin(M_PI*x)*cosh(M_PI*y);
        si = -cos(M_PI*x)*sinh(M_PI*y);
        q2 = sqrt(sr*sr+si*si);
        th2 = atan(si/sr);
        if (sr < 0.0) th2 += M_PI;
        gr = log(M_PI/(q1*q2))-gr;
        gi = -th1-th2-gi;
        x = x1;
        y = y1;
    }
    if (OPT == 0) {
        g0 = exp(gr);
        gr = g0*cos(gi);
        gi = g0*sin(gi);
    }
    g = std::complex<double>(gr,gi);
    return g;
}



void Rates::Initialize(tools::Property *options) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options, "xtp" );
    std::string key = "options." + Identify();

    /* ---- OPTIONS.XML Structure -----
     *
     * <rates>
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

    tools::vec F = options->get(key+".field").as<tools::vec>();
    _F = F;


    // Method
    if (options->exists(key+".method")) {
        _rateType = options->get(key+".method").as<std::string> ();
    }
    else {
        std::cout << std::endl 
             << "... ... ERROR: No method to calculate rates specified. "
             << std::endl;
        throw std::runtime_error("Missing input in options file.");
    }
    if (_rateType != "marcus" && _rateType != "jortner" && _rateType != "weissdorsey" && _rateType != "sven") {
        std::cout << std::endl
             << "... ... ERROR: Unknown rate type '" << _rateType << "' "
             << std::endl;
        throw std::runtime_error("Faulty input in options file.");
    }


    if (_rateType == "jortner") {
        
   
    // Vibrational quanta
    if (options->exists(key+".nmaxvib")) {
            _nMaxVib = options->get(key+".nmaxvib").as<double> ();
        }
        else {
            _nMaxVib = 20;
            std::cout << std::endl << "... ... WARNING: No cut-off number for QM vibrations "
                            "provided, using default 20.";
        }
        if (options->exists(key+".omegavib")) {
            _omegaVib = options->get(key+".omegavib").as<double> ();
        }
        else {
            _omegaVib = 0.2;
            std::cout << std::endl << "... ... WARNING: No QM vibration frequency provided, "
                            "using default 0.2eV.";
        }
    }

    
    if (_rateType == "weissdorsey") {
    // Kondo parameter
    if (options->exists(key+".kondo")) {
            _kondo = options->get(key+".kondo").as<double> ();
        }
        else {
            _kondo = 4.0;
            std::cout << std::endl << "... ... WARNING: No Kondo parameter provided. Using default 4.0.";
        }
    }

    // this->ParseEnergiesXML(top, options);

}


void Rates::ParseEnergiesXML(Topology *top, tools::Property *opt) {

    std::string key = "options.rates";
    std::string energiesXML = opt->get(key+".energiesXML").as<std::string> ();

    std::cout << std::endl
         << "... ... Site, reorg. energies from " << energiesXML << ". "
         << std::flush;

    tools::Property alloc;
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
    std::list<tools::Property*> mols = alloc.Select(key);
    std::list<tools::Property*> ::iterator molit;
    for (molit = mols.begin(); molit != mols.end(); ++molit) {

        key = "segments.segment";
        std::list<tools::Property*> segs = (*molit)->Select(key);
        std::list<tools::Property*> ::iterator segit;

        for (segit = segs.begin(); segit != segs.end(); ++segit) {

            std::string segName = (*segit)->get("name").as<std::string> ();

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


void Rates::EvaluatePair(Topology *top, QMPair *qmpair) {

    std::cout << "\r... ... Evaluating pair " << qmpair->getId()+1 << ". " << std::flush;

    bool pair_has_e = false;
    bool pair_has_h = false;
    bool pair_has_s = false;
    bool pair_has_t = false;

    std::string segName1 = qmpair->first->getName();
    std::string segName2 = qmpair->second->getName();

    pair_has_e = qmpair->isPathCarrier(-1);
    pair_has_h = qmpair->isPathCarrier(+1);
    pair_has_s = qmpair->isPathCarrier(+2);
    pair_has_t = qmpair->isPathCarrier(+3);

//    try {
//        pair_has_e = _seg_has_e.at(segName1) && _seg_has_e.at(segName2);
//        pair_has_h = _seg_has_h.at(segName1) && _seg_has_h.at(segName2);
//    }
//    catch (out_of_range) {
//        std::cout << std::endl << "... ... WARNING: No energy information for pair ["
//                     << segName1 << ", " << segName2 << "]. "
//                     << "Skipping... " << std::endl;
//
//        return;
//    }

    if (pair_has_e) {
        this->CalculateRate(top, qmpair, -1);
    }
    if (pair_has_h) {
        this->CalculateRate(top, qmpair, +1);
    }
    if (pair_has_s) {
        this->CalculateRate(top, qmpair, +2);
    }
    if (pair_has_t) {
        this->CalculateRate(top, qmpair, +3);
    }
}


void Rates::CalculateRate(Topology *top, QMPair *qmpair, int state) {

    const double NM2M    = 1.e-9;
    const double hbar_eV = 6.58211899e-16;

    Segment *seg1 = qmpair->first;
    Segment *seg2 = qmpair->second;

    double rate12 = 0.;                                       // 1->2

    double rate21 = 0.;                                       // 2->1

    double rate_symm12 = 0;
    double rate_symm21 = 0;
    //double measure = 0;
    double reorg12=0;
    double reorg21=0;
    double dG_Site=0;
    double dG_Field=0;
    
    if (state<2){
        reorg12  = seg1->getU_nC_nN(state)                 // 1->2
                        + seg2->getU_cN_cC(state);
        reorg21  = seg1->getU_cN_cC(state)                 // 2->1
                        + seg2->getU_nC_nN(state);
        dG_Site  = seg2->getU_cC_nN(state)                 // 1->2 == - 2->1
                        + seg2->getEMpoles(state)
                        - seg1->getU_cC_nN(state)
                        - seg1->getEMpoles(state);
        dG_Field = - state * _F * qmpair->R() * NM2M;
    }
    else if (state>=2){
        reorg12  = seg1->getU_nX_nN(state)                 // 1->2
                        + seg2->getU_xN_xX(state);
        reorg21  = seg1->getU_xN_xX(state)                 // 2->1
                        + seg2->getU_nX_nN(state);
        dG_Site  = seg2->getU_xX_nN(state)                 // 1->2 == - 2->1
                        + seg2->getEMpoles(state)
                        - seg1->getU_xX_nN(state)
                        - seg1->getEMpoles(state);       
    }
    
    
    
    
    
    double lOut     = qmpair->getLambdaO(state);              // 1->2 == + 2->1

    
          // 1->2 == - 2->1

    double J2 = qmpair->getJeff2(state);                      // 1->2 == + 2->1
    double dG = dG_Field + dG_Site;

    // +++++++++++++ //
    // JORTNER RATES //
    // +++++++++++++ //

    if (_rateType == "jortner" && lOut < 0.) {
        std::cout << std::endl
             << "... ... ERROR: Pair " << qmpair->getId() << " has negative "
                "outer-sphere reorganization energy. Cannot calculate Jortner "
                "rates. "
             << std::endl;
        throw std::runtime_error("");
    }
    else if (_rateType == "jortner" && lOut < 0.01) {
        std::cout << std::endl
             << "... ... WARNING: Pair " << qmpair->getId() << " has small "
                "outer-sphere reorganization energy (" << lOut << "eV). Could "
                "lead to over-estimated Jortner rates."
             << std::endl;
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
        
    // ++++++++++++ //
    // WEISS DORSEY RATES //
    // See Asadi et al. Nature Comm. 2708 (DOI: 10.1038/ncomms2708)
    // equation 6 
    // ++++++++++++ //

    } else if (_rateType == "weissdorsey") {
        
        _kondo = _kondo/2+1; // going from alpha to alpha'

        reorg12 = reorg12 + lOut;
        reorg21 = reorg21 + lOut;
        
        double characfreq12 = reorg12 /2 /_kondo/hbar_eV;
        double characfreq21 = reorg21 /2 /_kondo/hbar_eV;
        
        std::complex<double> M_I = std::complex<double>(0.0,1.0);
        
       rate12 = J2/pow(hbar_eV,2)/characfreq12
                * pow((hbar_eV*characfreq12/2/M_PI/_kT), (1-2*_kondo))
                * pow(std::abs(ccgamma(_kondo+M_I*(+dG/2/M_PI/_kT),1)),2)
                * pow(ccgamma(2*_kondo,0).real(), -1) * exp(+dG/2/_kT)
                * exp(-std::abs(dG)/hbar_eV/characfreq12); 

        rate21 = J2/pow(hbar_eV,2)/characfreq21
                * pow((hbar_eV*characfreq21/2/M_PI/_kT), (1-2*_kondo))
                * pow(std::abs(ccgamma(_kondo+M_I*(-dG/2/M_PI/_kT),1)),2)
                * pow(ccgamma(2*_kondo,0).real(), -1) * exp(-dG/2/_kT)
                * exp(-std::abs(dG)/hbar_eV/characfreq12);

        
        
    // ++++++++++++ //
    // SYMMETRIC RATES //
    // ++++++++++++ //
        
    } else if (_rateType == "sven") {
         
        reorg12 = reorg12 + lOut;
        reorg21 = reorg21 + lOut;
       
        rate12 = J2 / hbar_eV * sqrt( M_PI / (reorg12*_kT) )
                * exp( - (+dG + reorg12)*(+dG + reorg12) / (4*_kT*reorg12) );

        rate21 = J2 / hbar_eV * sqrt( M_PI / (reorg21*_kT) )
                * exp( - (-dG + reorg21)*(-dG + reorg21) / (4*_kT*reorg21) );
        
        double e1  = seg2->getU_cC_nN(state) + seg2->getEMpoles(state);    // 1->2 == - 2->1
        double e2  =  seg1->getU_cC_nN(state) + seg1->getEMpoles(state);       
        
//
        rate_symm12 = J2 / hbar_eV * sqrt( M_PI / (reorg12*_kT) )
                * exp( - ( dG*dG + reorg12*reorg12 - 2.*dG*reorg12 )  / (4*_kT*reorg12) ) ;

        rate_symm21 = J2 / hbar_eV * sqrt( M_PI / (reorg21*_kT) )
                * exp( - ( dG*dG + reorg21*reorg21 + 2.*dG*reorg21 )  / (4*_kT*reorg21) ) ;
        
        double _rate_symm12 = rate12 * exp( -e1 / _kT );
        double _rate_symm21 = rate21 * exp( -e2 / _kT );
        
        std::cout << " " << qmpair->Seg1()->getId() << " " << qmpair->Seg2()->getId() <<
                rate_symm12 << " " << " " << rate_symm21 << " " <<
                _rate_symm12 << " " <<  _rate_symm21 << std::endl;
    }

    qmpair->setRate12(rate12, state);
    qmpair->setRate21(rate21, state);
    qmpair->setIsPathCarrier(true, state);

}

int Rates::Factorial(int i) {
    int k= 1;
    for (int j=1; j<i; j++) { k = k*(j+1); }
    return k;
}

}}










#endif
