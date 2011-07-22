#ifndef _ESHUFFLE_H
#define	_ESHUFFLE_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>
#include <stdlib.h>
#include <math.h>
#include <list>

/** \brief Shuffle site energies to get rid of spatial correlations, while preserving the distribution width 

The site energies are randomly assigned to sites which keeps the distribution of site energies but destroys the spatial correlations present e.g. due to charge-dipole interactions. Therefore the new distribution of site-energy differences may differ from the old distribution of site-energy differences.

Callname: eshuffle
*/
class Eshuffle : public QMCalculator
{
public:
    Eshuffle() {};
    ~Eshuffle() {};

    const char *Description() { return "Shuffles site energies (Coulomb) to get rid of spatial correlations, while preserving the distribution width"; }

    bool EvaluateFrame(QMTopology *top);

};

inline bool Eshuffle::EvaluateFrame(QMTopology *top) {

    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    vector<QMCrgUnit *>::iterator itl;
    vector<double> energies;
    int i_en = 0;

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        energies.push_back( (*itl)->getDouble("energy_coulomb") );
    }
    
    random_shuffle ( energies.begin(), energies.end() );

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        (*itl)->setDouble("energy_coulomb", energies[i_en] );
        i_en++;
    }

    return true;
}

#endif	/* _ESHUFFLE_H */

