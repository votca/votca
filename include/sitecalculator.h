/* 
 * File:   sitecalculator.h
 * Author: vehoff
 *
 * Created on April 7, 2010, 4:50 PM
 */

#ifndef _SITECALCULATOR_H
#define	_SITECALCULATOR_H

#include "qmcalculator.h"

class SiteCalculator : public QMCalculator{
public:
    SiteCalculator();
    virtual ~SiteCalculator();

    bool EvaluateFrame(QMTopology *top);
    virtual void EvaluateSite(CrgUnit *crg) {};
};

inline bool SiteCalculator::EvaluateFrame(QMTopology *top){
    list<CrgUnit*>& crglist = top.crglist();
    for(list<CrgUnit*>::iterator iter = crglist.begin();iter!=crglist.end();++iter)
        EvaluateSite(*iter);
    return true;
}

#endif	/* _SITECALCULATOR_H */

