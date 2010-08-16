#include "polymerrates.h"
#include <votca/csg/interaction.h>

using namespace votca::csg;

void PolymerRates::Initialize(QMTopology *top, Property *options)
{   
    _nu = options->get("options.calc_rates.polymer.nu").as<double>();;
    
    _j[0]=_j[1] = 0.;
    if(options->exists("options.calc_rates.polymer.J0"))
        _j[0] = options->get("options.calc_rates.polymer.J0").as<double>();;
    if(options->exists("options.calc_rates.polymer.J1"))
        _j[1] = options->get("options.calc_rates.polymer.J1").as<double>();;
    
        _kT = options->get("options.calc_rates.thermal_energy").as<double>();
    _E = options->get("options.calc_rates.e_field").as<vec>();
}

bool PolymerRates::EvaluateFrame(QMTopology *top)
{
    InteractionContainer::iterator iter;

    for(    iter = top->BondedInteractions().begin();
            iter != top->BondedInteractions().end();
            ++iter)
    {
        Interaction *i = *iter;
        if(i->BeadCount() == 2) { // we found a bond
            QMBead *b1 = dynamic_cast<QMBead*>(top->getBead(i->getBeadId(0)));
            QMBead *b2 = dynamic_cast<QMBead*>(top->getBead(i->getBeadId(1)));
            CrgUnit *crg1 = b1->GetCrgUnit();
            CrgUnit *crg2 = b2->GetCrgUnit();

            if((crg1 == crg2)
                    || (crg1 == NULL )
                    || (crg2 == NULL ))
                continue;
            if(crg1->getId() > crg2->getId())
                swap(crg1,crg2);

            QMPair *pair = top->nblist().FindPair(crg1, crg2);
            if(!pair) {
                pair = new QMPair(crg1, crg2, top);
                top->nblist().AddPair(pair);
            }
            double cos_u= b1->U() * b2->U() / (abs(b1->U())*abs(b2->U()));
            double J = _j[0] + abs(_j[1]*cos_u);
            pair->setRate12(CalcRate(pair->first, pair->second, pair->r(), J));
            pair->setRate21(CalcRate(pair->second, pair->first, -pair->r(), J));
        }
    }
}

double PolymerRates::CalcRate(CrgUnit *crg1, CrgUnit *crg2, vec dist, double J)
{    
    double dE_nofield = crg2->getEnergy() - crg1->getEnergy();
    double reorg = 0.5*(crg1->getType()->getReorg() + crg2->getType()->getReorg());
    
    double dE =  dE_nofield - unit<nm,m>::to(dist) * _E;
    double DG_star = ( dE + reorg)*(dE + reorg)/(4*reorg);

    DG_star -= fabs(J);
    double rate = _nu * exp ( -DG_star  / _kT  ) ;
}
