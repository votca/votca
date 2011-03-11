/* 
 * File:   jortnerrates.h
 * Author: mayfalk
 *
 * Created on March 3, 2011, 4:12 PM
 */

#ifndef JORTNERRATES_H
#define	JORTNERRATES_H


#include "paircalculator.h"
#include <math.h>

class JortnerRates : public PairCalculator
{
public:
    JortnerRates() {};
    ~JortnerRates() {};

    void Initialize(QMTopology *top, Property *options);
    void EvaluatePair(QMTopology *top, QMPair *pair);

private:
    vec _E;
    double _kT;
    double _omegavib;
    int _nmaxvib;
     int Factorial(int i);
};

inline void JortnerRates::Initialize(QMTopology *top, Property *options){
    _kT = options->get("options.calc_rates.thermal_energy").as<double>();
    _E = options->get("options.calc_rates.e_field").as<vec>();
    _omegavib=0.2;
    _nmaxvib=5;
}

inline void JortnerRates::EvaluatePair(QMTopology *top, QMPair* pair){
    double rate_12 = 0.0;
    double rate_21 = 0.0;
    double Jeff2 = pair->calcJeff2();
    double lambda_outer = pair->getLambdaOuter();
    CrgUnit *crg1 = pair->first;
    CrgUnit *crg2 = pair->second;
    /// prefactor for future modifications
    double prefactor = 1.0;
    /// reorganization energy in eV as given in list_charges.xml
    double reorg = 0.5 * (crg1->getType()->getReorg()+crg2->getType()->getReorg());
    ///Huang Rhys
    double huang_rhys = reorg/_omegavib;
    /// free energy difference due to electric field, i.e. E*r_ij
    double dG_field = -_E * unit<nm,m>::to(pair->r());
    /// free energy difference due to different energy levels of molecules
    double dG_en = crg2->getEnergy() - crg1->getEnergy();
    /// electrostatics are taken into account in qmtopology and are contained in Energy
    /// total free energy difference
    double dG = dG_field + dG_en;
    
    int nn;
    for (nn = 0; nn<=_nmaxvib; nn++) {
        /// Jortner rate from first to second
    dG = dG_field + dG_en;
    rate_12 = rate_12 + prefactor * sqrt(M_PI/(lambda_outer * _kT)) * Jeff2 * exp(-huang_rhys) * pow(huang_rhys,nn) / Factorial(nn) *
            exp (-(dG + nn*_omegavib + lambda_outer)*(dG + nn*_omegavib +lambda_outer)/(4*_kT*lambda_outer))/hbar_eV;
    /// Jortner rate from second to first (dG_field -> -dG_field)
    dG = -dG_field - dG_en;
    rate_21 = rate_21 + prefactor * sqrt(M_PI/(lambda_outer * _kT)) * Jeff2 * exp(-huang_rhys) * pow(huang_rhys,nn) / Factorial(nn) *
            exp (-(dG + nn*_omegavib +lambda_outer)*(dG + nn*_omegavib +lambda_outer)/(4*_kT*lambda_outer))/hbar_eV;
    }
    pair->setRate12(rate_12);
    //cout<<rate_12;
    pair->setRate21(rate_21);
}

int JortnerRates::Factorial(int i) {
    int k,j;
k=1;
for(j=1;j<i;j++)
  {
  k=k*(j+1);
  }
return k;
}

#endif	/* JORTNERRATES_H */

