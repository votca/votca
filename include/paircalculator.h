/* 
 * File:   paircalculator.h
 * Author: vehoff
 *
 * Created on March 31, 2010, 5:58 PM
 */

#ifndef _PAIRCALCULATOR_H
#define	_PAIRCALCULATOR_H

#include "qmcalculator.h"

class PairCalculator : public QMCalculator{
public:
    PairCalculator();
    virtual ~PairCalculator();

    void EvaluateFrame(QMTopology *top);
    virtual void EvaluatePair(QMTopology *top, QMPair *pair) {};
};

void PairCalculator::EvaluateFrame(QMTopology *top){
    QMNBList &nblist = top.nblist();
    for(QMNBList::iterator iter = nblist.begin();iter!=nblist.end();++iter)
        EvaluatePair(*iter);
}

#endif	/* _PAIRCALCULATOR_H */

