/* 
 * File:   sitecalculator.h
 * Author: vehoff
 *
 * Created on April 7, 2010, 4:50 PM
 */

#ifndef _SITECALCULATOR_H
#define	_SITECALCULATOR_H

#include "qmcalculator.h"

namespace votca { namespace ctp {

class SiteCalculator : public QMCalculator{
public:
    SiteCalculator();
    virtual ~SiteCalculator();

    bool EvaluateFrame(QMTopology *top);
    virtual void EvaluateSite(QMCrgUnit *crg) {};
};

inline bool SiteCalculator::EvaluateFrame(QMTopology *top){
    vector<QMCrgUnit*>& crglist = top.CrgUnits();
    for(vector<QMCrgUnit*>::iterator iter = crglist.begin();iter!=crglist.end();++iter)
        EvaluateSite(*iter);
    return true;
}

}}

#endif	/* _SITECALCULATOR_H */

