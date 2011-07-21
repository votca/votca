#ifndef __VOTCA_MD2QM_OCCEQUILIBRIUM_H
#define	__VOTCA_MD2QM_OCCEQUILIBRIUM_H

#include <votca/ctp/qmcalculator.h>

class OccEquilibrium : public QMCalculator
{
public:
    OccEquilibrium() {}
    ~OccEquilibrium() {}

    const char *Description() { return "TODO"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
private:
    double _kT;
};

inline void OccEquilibrium::Initialize(QMTopology *top, Property *options)
{
    _kT = options->get("options.calc_rates.thermal_energy").as<double>();
}

inline bool OccEquilibrium::EvaluateFrame(QMTopology *top)
{
    double Ptot = 0;
    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    vector<QMCrgUnit *>::iterator itl;


    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        //double reorg = (*itl)->getType()->getReorg();
        double p = exp(-((*itl)->getTotalEnergy()) / _kT);
        Ptot +=p;
        (*itl)->setOccupationProbability(p);
    }
    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        (*itl)->setOccupationProbability(
                (*itl)->getOccupationProbability() / Ptot);
    }
}


#endif	/* OCCEQUILIBRIUM_H */

