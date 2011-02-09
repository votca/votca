/* 
 * File:   polymerrates.h
 * Author: ruehle
 *
 * Created on August 16, 2010, 12:13 PM
 */

#ifndef _POLYMERRATES_H
#define	_POLYMERRATES_H

#include "qmcalculator.h"

class PolymerRates : public QMCalculator
{
public:
    PolymerRates() {};
    ~PolymerRates() {};

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
//    void EndEvaluate(QMTopology *top);

    double CalcRate(QMCrgUnit *crg1, QMCrgUnit *crg2, vec dist, double J);

private:
    double _nu;
    double _J;
    double _kT;
    vec _E;
};

#endif	/* POLYMERRATES_H */

