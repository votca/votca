#ifndef _CALC_INTEGRALS_H
#define	_CALC_INTEGRALS_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/paircalculator.h>

/**
	\brief Semi-empirical electronic coupling elements for all neighbor list pairs

Semi-emprirical (ZINDO) electronic coupling elements for all conjuageted segments from the neighbout list (J. Int. J. Quantum Chem. 2008, 108, 51-56.). Requires molecular orbitals in GAUSSIAN format.

Callname: izindo    

*/

class Izindo : public PairCalculator
{
public:
    Izindo() {};
    ~Izindo() {};

    const char *Description() { return "Semi-empirical electronic coupling elements for all neighbor list pairs"; }

    void EvaluatePair(QMTopology *top, QMPair *pair);
};

inline void Izindo::EvaluatePair(QMTopology *top, QMPair *pair){
    CrgUnit *crg1 = pair->Crg1PBCCopy();
    CrgUnit *crg2 = pair->Crg2PBCCopy();
    vector <double> Js = top->GetJCalc().CalcJ(*crg1, *crg2);
    pair->setJs(Js);
}

#endif	/* _CALC_INTEGRALS_H */

