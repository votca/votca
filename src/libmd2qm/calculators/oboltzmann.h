#ifndef OBOLTZMANN_H
#define OBOLTZMANN_H

#include <votca/ctp/qmcalculator.h>
/**
        \brief Site occupations as a Boltzmann distribution of site energies.

For an ergodic system in equilibrium the site occupation can be calculated from the site energy using the Boltzmann distribution
\f[ p_i = \exp \left(- \frac{E_i}{k_\textrm{B} T} \right) / \sum_i  \exp \left(- \frac{E_i}{k_\textrm{B} T} \right) \f].

Callname: oboltzmann

*/
class Oboltzmann : public QMCalculator
{
public:
    Oboltzmann() {}
    ~Oboltzmann() {}

    const char *Description() { return "Site occupations as a Boltzmann distribution of site energies"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
private:
    double _kT;
};

inline void Oboltzmann::Initialize(QMTopology *top, Property *options)
{
     double temp= options->get("options.calc_rates.temperature").as<double>();
    _kT=temp*8.6173324e-5;
}

inline bool Oboltzmann::EvaluateFrame(QMTopology *top)
{
    double Ptot = 0;
    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    vector<QMCrgUnit *>::iterator itl;


    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        double p = exp(-((*itl)->getTotalEnergy()) / _kT);
        Ptot +=p;
        (*itl)->setOccupationProbability(p);
    }
    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        (*itl)->setOccupationProbability(
                (*itl)->getOccupationProbability() / Ptot);
    }
}


#endif	/* OBOLTZMANN_H */

