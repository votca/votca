/* 
 * File:   calc_integrals.h
 * Author: vehoff
 *
 * Created on April 8, 2010, 10:33 AM
 */

#ifndef _CALC_INTEGRALS_H
#define	_CALC_INTEGRALS_H

#include "qmpair.h"
#include "paircalculator.h"

class CalcIntegrals : public PairCalculator
{
public:
    CalcIntegrals() {};
    ~CalcIntegrals() {};

    void EvaluatePair(QMTopology *top, QMPair *pair);
};

void CalcIntegrals::EvaluatePair(QMTopology *top, QMPair *pair){
    CrgUnit *crg1 = pair->Crg1();
    CrgUnit *crg2 = pair->Crg2();
    vector <double> Js = top->GetJCalc().CalcJ(*crg1, *crg2);
    pair->setJs(Js);
}

#endif	/* _CALC_INTEGRALS_H */

