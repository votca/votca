#ifndef _CALC_INTEGRALS_H
#define	_CALC_INTEGRALS_H

#include "qmpair.h"
#include "paircalculator.h"

/**
	\brief Semi-empirical electronic coupling elements for all neighbor list pairs

Semi-emprirical (ZINDO) electronic coupling elements for all conjuageted segments from the neighbout list. Requires molecular orbitals in GAUSSIAN format.

Callname: izindo    

References: Kirkpatrick, J. Int. J. Quantum Chem. 2008, 108, 51-56.
*/

class CalcIntegrals : public PairCalculator
{
public:
    CalcIntegrals() {};
    ~CalcIntegrals() {};

    const char *Description() { return "Semi-empirical electronic coupling elements for all neighbor list pairs"; }

    void EvaluatePair(QMTopology *top, QMPair *pair);
};

inline void CalcIntegrals::EvaluatePair(QMTopology *top, QMPair *pair){
    CrgUnit *crg1 = pair->Crg1();
    CrgUnit *crg2 = pair->Crg2();
    vector <double> Js = top->GetJCalc().CalcJ(*crg1, *crg2);
    pair->setJs(Js);
}

#endif	/* _CALC_INTEGRALS_H */

