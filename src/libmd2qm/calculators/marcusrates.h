#ifndef _MARCUS_RATES_H
#define	_MARCUS_RATES_H

#include "paircalculator.h"

/**
    \brief Calculates hopping rates using Marcus theroy

Callname: marcusrates

Hopping rates between all pairs i,j are computed from the high temperature limit of Marcus theory
    \f[\omega_{ij}  = \frac{2 \pi}{\hbar}  \frac{ J_{ij}^2 }{\sqrt{ 4 \pi \lambda_{ij} k_B T}} \exp \left[-\frac{\left(\Delta E_{ij}-\lambda_{ij}\right)^2}{4 \lambda_{ij}k_B T} \right]
    \f]
where \f$T\f$ is the temperature, \f$\lambda_{ij} = \lambda_{ij}^{int} + \lambda_{ij}^{out}\f$ is the reorganization energy, which is a sum of intra- and inter-molecular (outer-sphere) contributions, \f$\Delta E_{ij}\f$ is the site-energy difference, or driving force, and \f$J_{ij}\f$ is the electronic coupling element, or transfer integral. You should have computed transfer integrals and site energies before.   

Reference: Marcus, R. A. Reviews of Modern Physics 1993, 65, 599. 
*/

class MarcusRates : public PairCalculator
{
public:
    MarcusRates() {};
    ~MarcusRates() {};

    void Initialize(QMTopology *top, Property *options);
    void EvaluatePair(QMTopology *top, QMPair *pair);

private:
    vec _E;
    double _kT;
};

inline void MarcusRates::Initialize(QMTopology *top, Property *options){
    _kT = options->get("options.calc_rates.thermal_energy").as<double>();
    _E = options->get("options.calc_rates.e_field").as<vec>();
}

inline void MarcusRates::EvaluatePair(QMTopology *top, QMPair* pair){
    double rate_12 = 0.0;
    double rate_21 = 0.0;
    double Jeff2 = pair->calcJeff2();
    CrgUnit *crg1 = pair->first;
    CrgUnit *crg2 = pair->second;
    /// prefactor for future modifications
    double prefactor = 1.0;
    /// reorganization energy in eV as given in list_charges.xml
    double reorg = 0.5 * (crg1->getType()->getReorg()+crg2->getType()->getReorg());
    /// free energy difference due to electric field, i.e. E*r_ij
    double dG_field = -_E * unit<nm,m>::to(pair->r());
    /// free energy difference due to different energy levels of molecules
    double dG_en = crg2->getEnergy() - crg1->getEnergy();
    /// electrostatics are taken into account in qmtopology and are contained in Energy
    /// total free energy difference
    double dG = dG_field + dG_en;
    /// Marcus rate from first to second
    rate_12 = prefactor * sqrt(M_PI/(reorg * _kT)) * Jeff2 *
            exp (-(dG + reorg)*(dG + reorg)/(4*_kT*reorg))/hbar_eV;
    /// Marcus rate from second to first (dG_field -> -dG_field)
    dG = -dG_field - dG_en;
    rate_21 = prefactor * sqrt(M_PI/(reorg * _kT)) * Jeff2 *
            exp (-(dG + reorg)*(dG + reorg)/(4*_kT*reorg))/hbar_eV;

    pair->setRate12(rate_12);
    pair->setRate21(rate_21);
}

#endif	/* _MARCUS_RATES_H */

