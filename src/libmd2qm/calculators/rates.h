#ifndef RATES_H
#define	RATES_H


#include <votca/ctp/paircalculator.h>
#include <math.h>

/**
    \brief Calculates hopping rates using either Marcu's or Jortner's derivation

Callname: rates

Hopping rates between all pairs i,j are computed.
This can be done using the high temperature limit of Marcus theory
    \f[\omega_{ij}  = \frac{2 \pi}{\hbar}  \frac{ J_{ij}^2 }{\sqrt{ 4 \pi \lambda_{ij} k_B T}} \exp \left[-\frac{\left(\Delta E_{ij}-\lambda_{ij}\right)^2}{4 \lambda_{ij}k_B T} \right]
    \f]
where \f$T\f$ is the temperature, \f$\lambda_{ij} = \lambda_{ij}^{int} + \lambda_{ij}^{out}\f$ is the reorganization energy, which is a sum of intra- and inter-molecular (outer-sphere) contributions, \f$\Delta E_{ij}\f$ is the site-energy difference, or driving force, and \f$J_{ij}\f$ is the electronic coupling element, or transfer integral. You should have computed transfer integrals and site energies before.   

Alternatively, Jortner's approach can used by treating both molecules quantum mechanically in presence of a classical outer sphere mode 
    \f[\omega_{ij}  = \frac{2 \pi}{\hbar}  \frac{ J_{ij}^2 }{\sqrt{ 4 \pi \lambda_{ij}^{out} k_B T}} \sum_{N=0}^\infty \frac{1}{N!} \left( \frac{\lambda_{ij}^{int}}{\hbar\omega^{int}} \right)^{N} \exp \left( - \frac{\lambda_{ij}^{int}}{\hbar\omega^{int}}\right) 
\exp \left\{ -\frac{ \left[ \Delta E_{ij}-\hbar N\omega^{int} -\lambda_{ij}^{out} \right]^2}{4\lambda_{ij}^{out} k_B T}\right\}\f]
where \f$T\f$ is the temperature, \f$\lambda_{ij}^{int} , \lambda_{ij}^{out}\f$ are the intramolecular and outer sphere reorganization energies,\f$\omega^{int}\f$ is the averaged vibrational frequency of the intramolecular modes (usuallay assumed to be close to the C-C bond stretch at 0.2eV). \f$\Delta E_{ij}\f$ is the site-energy difference, or driving force, and \f$J_{ij}\f$ is the electronic coupling element, or transfer integral. You should have computed transfer integrals, site energies and outer sphere reorganization energies before.  

Reference: V. May and O. Kuhn, Charge and Energy Transfer Dynamics in Molecular Systems. 
*/


class Rates : public PairCalculator
{
public:
    Rates() {};
    ~Rates() {};

    const char *Description() { return "Calculates hopping rates using Jortner's derivation"; }

    void Initialize(QMTopology *top, Property *options);
    void EvaluatePair(QMTopology *top, QMPair *pair);

private:
    char _rate_type; //M=Marcus J=Jortner
    vec _E;
    double _kT;
    double _omegavib;
    int _nmaxvib;
     int Factorial(int i);
};

inline void Rates::Initialize(QMTopology *top, Property *options){
	
    _kT = options->get("options.calc_rates.thermal_energy").as<double>();
    _E = options->get("options.calc_rates.e_field").as<vec>();
   // _omegavib=0.2;
   // _nmaxvib=10;
      if (options->exists("options.calc_rates.method")) {
        if (options->get("options.calc_rates.method").as<string > () == "marcus") {
		_rate_type='M';
            cout << "Computing rates from Marcus theory"<< endl;
        }
        else if (options->get("options.calc_rates.method").as<string > () == "jortner") {
		_rate_type='J';
        if (options->exists("options.calc_rates.nmaxvib")) {
        _nmaxvib = options->get("options.calc_rates.nmaxvib").as<double>();
    } else {
        _nmaxvib = 20;
        cout << "Warning: no cutoff number for qm vibrations  provided, using default 20" << endl;
    }
        if (options->exists("options.calc_rates.omegavib")) {
        _omegavib = options->get("options.calc_rates.omegavib").as<double>();
    } else {
        _omegavib = 0.2;
        cout << "Warning: no qm vibration frequency  provided, using default 0.2eV" << endl;
    }

            cout << "Computing rates from Jortner summing up to " <<_nmaxvib << " quanta using a frequency of "<<_omegavib << " eV" << endl;
        }

        else throw std::runtime_error("Error in CalcRates::Initialize : no such rates method, should be marcus  or jortner");
    } else throw std::runtime_error("Error in Rates:Initialize : no rates_method specified");

}

inline void Rates::EvaluatePair(QMTopology *top, QMPair* pair){
    double rate_12 = 0.0;
    double rate_21 = 0.0;
    double Jeff2 = pair->calcJeff2();
    if (pair->DoubleExists("lambda_outer")){
    double lambda_outer = pair->getDouble("lambda_outer");
}else {if (_rate_type=='J')  throw std::runtime_error("Error in CalcRates::EvaluatePair : this pair has not outer sphere reorganization energy necessary to compute Jortner rates. Compute lambdaouter or use Marcus rates.");
}

    QMCrgUnit *crg1 = pair->first;
    QMCrgUnit *crg2 = pair->second;
    /// prefactor for future modifications
    double prefactor = 1.0;
    /// reorganization energy in eV as given in list_charges.xml
    double reorg12 = crg1->getDouble("lambda_intra_discharging")+crg2->getDouble("lambda_intra_charging");
    double reorg21 = crg2->getDouble("lambda_intra_discharging")+crg1->getDouble("lambda_intra_charging");
    /// free energy difference due to electric field, i.e. E*r_ij
    double dG_field = -_E * unit<nm,m>::to(pair->r());
    /// free energy difference due to different energy levels of molecules
    double dG_en = crg2->getTotalEnergy() - crg1->getTotalEnergy();
    /// electrostatics are taken into account in qmtopology and are contained in Energy
    /// total free energy difference
    double dG = dG_field + dG_en;
 
       if (_rate_type=='J'){ 
    ///Huang Rhys
    double huang_rhys12 = reorg12/_omegavib;
    double huang_rhys21 = reorg21/_omegavib;
 
    int nn;
    for (nn = 0; nn<=_nmaxvib; nn++) {
        /// Jortner rate from first to second
    dG = dG_field + dG_en;
    rate_12 = rate_12 + prefactor * sqrt(M_PI/(pair->getDouble("lambda_outer") * _kT)) * Jeff2 * exp(-huang_rhys12) * pow(huang_rhys12,nn) / Factorial(nn) *
            exp (-(dG + nn*_omegavib + pair->getDouble("lambda_outer"))*(dG + nn*_omegavib +pair->getDouble("lambda_outer"))/(4*_kT*pair->getDouble("lambda_outer")))/hbar_eV;
    /// Jortner rate from second to first (dG_field -> -dG_field)
    dG = -dG_field - dG_en;
    rate_21 = rate_21 + prefactor * sqrt(M_PI/(pair->getDouble("lambda_outer") * _kT)) * Jeff2 * exp(-huang_rhys21) * pow(huang_rhys21,nn) / Factorial(nn) *
            exp (-(dG + nn*_omegavib +pair->getDouble("lambda_outer"))*(dG + nn*_omegavib +pair->getDouble("lambda_outer"))/(4*_kT*pair->getDouble("lambda_outer")))/hbar_eV;
    }
}//end J

if (_rate_type=='M'){
 if (pair->DoubleExists("lambda_outer")){
reorg12=reorg12+ pair->getDouble("lambda_outer");
reorg21=reorg21+ pair->getDouble("lambda_outer");}
rate_12 = prefactor * sqrt(M_PI/(reorg12 * _kT)) * Jeff2 *
            exp (-(dG + reorg12)*(dG + reorg12)/(4*_kT*reorg12))/hbar_eV;
    /// Marcus rate from second to first (dG_field -> -dG_field)
    dG = -dG_field - dG_en;
    rate_21 = prefactor * sqrt(M_PI/(reorg21 * _kT)) * Jeff2 *
            exp (-(dG + reorg21)*(dG + reorg21)/(4*_kT*reorg21))/hbar_eV;
}//end M
    pair->setRate12(rate_12);
    //cout<<rate_12;
    pair->setRate21(rate_21);
}

int Rates::Factorial(int i) {
    int k,j;
k=1;
for(j=1;j<i;j++)
  {
  k=k*(j+1);
  }
return k;
}

#endif	/* RATES_H */

