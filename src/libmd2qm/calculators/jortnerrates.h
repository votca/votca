#ifndef JORTNERRATES_H
#define	JORTNERRATES_H


#include <votca/ctp/paircalculator.h>
#include <math.h>

/**
    \brief Calculates hopping rates using Jortner's derivation

Callname: jortnerrates

Hopping rates between all pairs i,j are computed by treating both molecules quantum mechanically in presence of a classical outer sphere mode 
    \f[\omega_{ij}  = \frac{2 \pi}{\hbar}  \frac{ J_{ij}^2 }{\sqrt{ 4 \pi \lambda_{ij}^{out} k_B T}} \sum_{N=0}^\infty \frac{1}{N!} \left( \frac{\lambda_{ij}^{int}}{\hbar\omega^{int}} \right)^{N} \exp \left( - \frac{\lambda_{ij}^{int}}{\hbar\omega^{int}}\right) 
\exp \left\{ -\frac{ \left[ \Delta E_{ij}-\hbar N\omega^{int} -\lambda_{ij}^{out} \right]^2}{4\lambda_{ij}^{out} k_B T}\right\}\f]
where \f$T\f$ is the temperature, \f$\lambda_{ij}^{int} , \lambda_{ij}^{out}\f$ are the intramolecular and outer sphere reorganization energies,\f$\omega^{int}\f$ is the averaged vibrational frequency of the intramolecular modes (usuallay assumed to be close to the C-C bond stretch at 0.2eV). \f$\Delta E_{ij}\f$ is the site-energy difference, or driving force, and \f$J_{ij}\f$ is the electronic coupling element, or transfer integral. You should have computed transfer integrals, site energies and outer sphere reorganization energies before.  

Reference: V. May and O. Kuhn, Charge and Energy Transfer Dynamics in Molecular Systems. 
*/


class JortnerRates : public PairCalculator
{
public:
    JortnerRates() {};
    ~JortnerRates() {};

    const char *Description() { return "Calculates hopping rates using Jortner's derivation"; }

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
    _nmaxvib=10;
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

