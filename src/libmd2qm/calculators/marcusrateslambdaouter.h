/* 
 * File:   MarcusRatesLambdaOutlambdaout.h
 * Author: mayfalk
 *
 * Created on March 9, 2011, 10:59 AM
 */

/*
 * File:   marcus_rates.h
 * Author: vehoff
 *
 * Created on April 8, 2010, 10:53 AM
 */

#ifndef _MARCUS_RATES_LAMBDA_OUTER_H
#define	_MARCUS_RATES_LAMBDA_OUTER_H

#include <votca/ctp/paircalculator.h>

class MarcusRatesLambdaOuter : public PairCalculator
{
public:
    MarcusRatesLambdaOuter() {};
    ~MarcusRatesLambdaOuter() {};

    const char *Description() { return "TODO: give description"; }

    void Initialize(QMTopology *top, Property *options);
    void EvaluatePair(QMTopology *top, QMPair *pair);

private:
    vec _E;
    double _kT;
};

inline void MarcusRatesLambdaOuter::Initialize(QMTopology *top, Property *options){
    _kT = options->get("options.calc_rates.thermal_energy").as<double>();
    _E = options->get("options.calc_rates.e_field").as<vec>();
}

inline void MarcusRatesLambdaOuter::EvaluatePair(QMTopology *top, QMPair* pair){
    double rate_12 = 0.0;
    double rate_21 = 0.0;
    double Jeff2 = pair->calcJeff2();
    QMCrgUnit *crg1 = pair->first;
    QMCrgUnit *crg2 = pair->second;
    /// prefactor for future modifications
    double prefactor = 1.0;
    /// reorganization energy in eV as given in list_charges.xml
    double reorg12 = crg1->getDouble("lambda_intra_discharging")+crg2->getDouble("lambda_intra_charging");
    double reorg21 = crg2->getDouble("lambda_intra_discharging")+crg1->getDouble("lambda_intra_charging");
    /// outer sphere lambda
    double lambda_outer = pair->getDouble("lambda_outer");
    ///
    reorg21=reorg21+lambda_outer;
    reorg12=reorg12+lambda_outer;
    /// free energy difference due to electric field, i.e. E*r_ij
    double dG_field = -_E * unit<nm,m>::to(pair->r());
    /// free energy difference due to different energy levels of molecules
    double dG_en = crg2->getTotalEnergy() - crg1->getTotalEnergy();
    /// electrostatics are taken into account in qmtopology and are contained in Energy
    /// total free energy difference
    double dG = dG_field + dG_en;
    /// Marcus rate from first to second
    rate_12 = prefactor * sqrt(M_PI/(reorg12 * _kT)) * Jeff2 *
            exp (-(dG + reorg12)*(dG + reorg12)/(4*_kT*reorg12))/hbar_eV;
    /// Marcus rate from second to first (dG_field -> -dG_field)
    dG = -dG_field - dG_en;
    rate_21 = prefactor * sqrt(M_PI/(reorg21 * _kT)) * Jeff2 *
            exp (-(dG + reorg21)*(dG + reorg21)/(4*_kT*reorg21))/hbar_eV;

    pair->setRate12(rate_12);
    pair->setRate21(rate_21);
}

#endif	/* _MARCUS_RATES_LAMBDA_OUTER_H */

