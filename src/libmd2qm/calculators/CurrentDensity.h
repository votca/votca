/* 
 * File:   CurrentDensity.h
 * Author: lukyanov
 *
 * Created on March 2, 2011, 1:06 PM
 */

#ifndef CURRENTDENSITY_H
#define	CURRENTDENSITY_H

#include "qmpair.h"
#include "paircalculator.h"

class CurrentDensity : public PairCalculator
{
public:
    CurrentDensity() {};
    ~CurrentDensity() {};

    void EvaluatePair(QMTopology *top, QMPair *pair);
    void EndEvaluate(QMTopology *top);

private:
    map <QMCrgUnit *, vec> _map_crgunt_current;
};

inline void CurrentDensity::EvaluatePair(QMTopology *top, QMPair *pair){
    QMCrgUnit *crg1 = pair->Crg1();
    QMCrgUnit *crg2 = pair->Crg2();

    double p1 = crg1->getOccupationProbability();
    double p2 = crg2->getOccupationProbability();

    double w12 = pair->rate12();
    double w21 = pair->rate21();
    vec r = pair->r();

    _map_crgunt_current[crg1] += 0.5 * (p2 * w21 - p1 * w12) * r;
    _map_crgunt_current[crg2] -= 0.5 * (p2 * w21 - p1 * w12) * r;

}


void CurrentDensity::EndEvaluate(QMTopology* top) {
    // normalize
    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    vector<QMCrgUnit *>::iterator itl;
    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
      cout <<  (*itl)->getId() << "\t" << (*itl)->GetCom()<< "\t" << _map_crgunt_current[(*itl)] << endl;
    }
}

#endif	/* CURRENTDENSITY_H */

