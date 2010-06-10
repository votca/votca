/* 
 * File:   Estatic_Calculator.h
 * Author: mayfalk
 *
 * Created on May 21, 2010, 11:23 AM
 */

#ifndef _CALC_ESTATICS_H
#define	_CALC_ESTATICS_H

#include "qmpair.h"
#include "qmcalculator.h"


class CalcEstatics : public QMCalculator
{
public:
    CalcEstatics() {};
    ~CalcEstatics() {};

    double CalcPot(Topology *atop, Molecule *mol);
    double CalcPot2(Topology *atop, Molecule *mol);
    void EvaluateFrame(QMTopology *top);
};

#endif	/* _CALC_ESTATICS_H */

