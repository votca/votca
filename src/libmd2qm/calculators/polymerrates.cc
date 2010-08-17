#include "polymerrates.h"
#include <votca/csg/interaction.h>

using namespace votca::csg;

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

/*    for(    iter = top->BondedInteractions().begin();
            iter != top->BondedInteractions().end();
            ++iter)
    {
        Interaction *i = *iter;
        if(i->BeadCount() == 2) { // we found a bond
             QMBead *b1 = dynamic_cast<QMBead*>(top->getBead(i->getBeadId(0)));
            QMBead *b2 = dynamic_cast<QMBead*>(top->getBead(i->getBeadId(1)));

 */
    // iterator over all molecules
    for(MoleculeContainer::iterator imol = top->Molecules().begin();
          imol!=top->Molecules().end(); ++imol) {
        // now go over all connections and check weather they connect two charge units
        for(int i=0; i<(*imol)->BeadCount()-1; ++i) { // hack for linear molecules

            QMBead *b1 = dynamic_cast<QMBead*>((*imol)->getBead(i));
            QMBead *b2 = dynamic_cast<QMBead*>((*imol)->getBead(i+1));

            CrgUnit *crg1 = b1->GetCrgUnit();
            CrgUnit *crg2 = b2->GetCrgUnit();

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
double PolymerRates::CalcRate(CrgUnit *crg1, CrgUnit *crg2, vec dist, double J)
{    
    double dE_nofield = crg2->getEnergy() - crg1->getEnergy();
    double reorg = 0.5*(crg1->getType()->getReorg() + crg2->getType()->getReorg());
    
    double dE =  dE_nofield - unit<nm,m>::to(dist) * _E;
    double DG_star = ( dE + reorg)*(dE + reorg)/(4*reorg);

    // J lowers the barrier
    DG_star -= fabs(J);
    double rate = _nu * exp ( -DG_star  / _kT  ) ;
}
