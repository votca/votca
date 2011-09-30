/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef RATES_H
#define	RATES_H


#include <votca/ctp/paircalculator.h>
#include <math.h>

namespace votca { namespace ctp {
/**
    \brief Calculates hopping rates using derivations from Marcus or Jortner

Callname: rates

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

/* polymer rates [development]
private:
    double _nu;
    double _J;
    double _kT;
    vec _E;
*/

};

inline void Rates::Initialize(QMTopology *top, Property *options){
	
    double temp= options->get("options.calc_rates.temperature").as<double>();
    _kT=temp*8.6173324e-5;
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
    QMCrgUnit *crg1 = pair->first;
    QMCrgUnit *crg2 = pair->second;
            if (pair->DoubleExists("lambda_outer")) {
                double lambda_outer = pair->getDouble("lambda_outer");
                if ((_rate_type == 'J') && (lambda_outer < 0.0)) {
                    cout << "Pair" << crg1->getId() << " : " << crg2->getId() << " has negative outer sphere reorganization energy preventing computation  of Jortner rates." << endl;
                    throw std::runtime_error("Error in CalcRates::EvaluatePair: negative outer sphere reorganization");
                }
                if ((_rate_type == 'J') && (lambda_outer < 0.01)) 
                    cout << "Warning: in CalcRates::EvaluatePair : pair " << crg1->getId() << " : " << crg2->getId() << " has very small outer sphere reorganization energy (< 0.01eV) which might produce too high Jortner rates."<<endl;
           
            }
            else {
                if (_rate_type == 'J') {
                    cout << "Pair" << crg1->getId() << " : " << crg2->getId() << ": has not outer sphere reorganization energy necessary to compute Jortner rates. Compute lambdaouter or use Marcus rates."<<endl;
                            throw std::runtime_error("Error in CalcRates::EvaluatePair: no outer sphere reorganization");
                }
            }

    
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


/* polymer rates (adiabatic rates from polypyrrole paper)

void PolymerRates::Initialize(QMTopology *top, Property *options)
{   
    // prefactor (something like a nuclear frequency)
    _nu = options->get("options.calc_rates.polymer.nu").as<double>();;
    
    _J=0;
    
    // transfer integrals for J=J_0 * cos(phi)
    if(options->exists("options.calc_rates.polymer.J0"))
        _J = options->get("options.calc_rates.polymer.J0").as<double>();;
    
    // thermal energy in eV
    _kT = options->get("options.calc_rates.thermal_energy").as<double>();
    // electric field vector in V/m
    _E = options->get("options.calc_rates.e_field").as<vec>();
}

bool PolymerRates::EvaluateFrame(QMTopology *top)
{
    InteractionContainer::iterator iter;

    // iterator over all molecules
    for(MoleculeContainer::iterator imol = top->Molecules().begin();
          imol!=top->Molecules().end(); ++imol) {
        // now go over all connections and check weather they connect two charge units
        for(int i=0; i<(*imol)->BeadCount()-1; ++i) { // hack for linear molecules

            QMBead *b1 = dynamic_cast<QMBead*>((*imol)->getBead(i));
            QMBead *b2 = dynamic_cast<QMBead*>((*imol)->getBead(i+1));

            QMCrgUnit *crg1 = b1->GetCrgUnit();
            QMCrgUnit *crg2 = b2->GetCrgUnit();

            // are the two beads on same crg unit or crgunit is null -> ignore
            if((crg1 == crg2)
                    || (crg1 == NULL )
                    || (crg2 == NULL ))
                continue;

            // swap if ordering is mixed up
            if(crg1->getId() > crg2->getId())
                swap(crg1,crg2);

            // is there already a pair for these two crgunits
            QMPair *pair = top->nblist().FindPair(crg1, crg2);
            if(!pair) { // if not create one
                pair = new QMPair(crg1, crg2, top);
                top->nblist().AddPair(pair);
            }
            
            // calculate the transfer integral
            double cos_u= b1->U() * b2->U() / (abs(b1->U())*abs(b2->U()));
            double J = _J*cos_u;

            // calculate rates for both directions
            pair->setRate12(CalcRate(pair->first, pair->second, pair->r(), J));
            pair->setRate21(CalcRate(pair->second, pair->first, -pair->r(), J));
        }
    }
}

// adiabatic rate expression
double PolymerRates::CalcRate(QMCrgUnit *crg1, QMCrgUnit *crg2, vec dist, double J)
{    
    double dE_nofield = crg2->getTotalEnergy() - crg1->getTotalEnergy();
    double reorg = 0.5*(crg1->getDouble("lambda_intra_discharging") + crg2->getDouble("lambda_charging")+crg2->getDouble("lambda_intra_discharging") + crg1->getDouble("lambda_charging"));
    
    double dE =  dE_nofield - unit<nm,m>::to(dist) * _E;
    double DG_star = ( dE + reorg)*(dE + reorg)/(4*reorg);

    // J lowers the barrier
    DG_star -= fabs(J);
    double rate = _nu * exp ( -DG_star  / _kT  ) ;
}
*/
}}

#endif	/* RATES_H */

