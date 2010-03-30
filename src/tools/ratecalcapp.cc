#include "ratecalcapp.h"

RateCalculator::RateCalculator()
{}

RateCalculator::~RateCalculator()
{}

void RateCalculator::HelpText(){
    cout << "CTP Ratecalculator." << endl;
}

void RateCalculator::Initialize(){
    _kT = _options.get("options.calc_rates.thermal_energy").as<double>();
    _E = _options.get("options.calc_rates.e_field").as<vec>();
}

bool RateCalculator::EvaluateFrame(){
    QMNBList &nblist = _qmtop.nblist();
    for(QMNBList::iterator iter = nblist.begin();iter!=nblist.end();++iter)
    {
        double rate_12 = 0.0;
        double rate_21 = 0.0;
        double Jeff2 = (*iter)->calcJeff2();
        CrgUnit *crg1 = (*iter)->first;
        CrgUnit *crg2 = (*iter)->second;
        /// prefactor for future modifications
        double prefactor = 1.0;
        /// reorganization energy in eV as given in list_charges.xml
        double reorg = 0.5 * (crg1->getType()->getReorg()+crg2->getType()->getReorg());
        /// free energy difference due to electric field, i.e. E*r_ij
        double dG_field = -_E * unit<nm,m>::to((*iter)->r());
        /// free energy difference due to different energy levels of molecules
        double dG_en = crg2->getEnergy() - crg1->getEnergy();
        /// electrostatics are taken into account in qmtopology and are contained in Energy
        /// total free energy difference
        double dG = dG_field + dG_en;
        /// Marcus rate from first to second
        rate_12 = prefactor * sqrt(M_PI/(reorg * _kT)) * Jeff2 *
                exp (-(dG + reorg)*(dG + reorg)/(4*_kT*reorg))/hbar_eV;
        /// Marcus rate from second to first (dG_field -> -dG_field)
        dG = -dG_field + dG_en;
        rate_21 = prefactor * sqrt(M_PI/(reorg * _kT)) * Jeff2 *
                exp (-(dG + reorg)*(dG + reorg)/(4*_kT*reorg))/hbar_eV;
        //cout << "rate_21 = " << rate_21 << endl;
        
        (*iter)->setRate12(rate_12);
        (*iter)->setRate21(rate_21);
    }
}

void RateCalculator::EndEvaluate()
 {
    stringstream ss;
    ss << _qmtop.getStep();
    string res;
    ss >> res;
    res = string ("nbl_") + res + string(".res");
    PrintNbs(res);
 }