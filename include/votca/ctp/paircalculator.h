/* 
 * File:   paircalculator.h
 * Author: vehoff
 *
 * Created on March 31, 2010, 5:58 PM
 */

#ifndef _PAIRCALCULATOR_H
#define	_PAIRCALCULATOR_H

#include "qmcalculator.h"

namespace votca { namespace ctp {

class PairCalculator : public QMCalculator{
public:
    PairCalculator() {};
    virtual ~PairCalculator() {};

    bool EvaluateFrame(QMTopology *top);
    virtual void EvaluatePair(QMTopology *top, QMPair *pair) {};
};

inline bool PairCalculator::EvaluateFrame(QMTopology *top){
    QMNBList &nblist = top->nblist();
    for(QMNBList::iterator iter = nblist.begin();iter!=nblist.end();++iter)
        EvaluatePair(top, *iter);
    return true;
}

}}

#endif	/* _PAIRCALCULATOR_H */

