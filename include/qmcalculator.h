/* 
 * File:   qmcalculator.h
 * Author: vehoff
 *
 * Created on March 31, 2010, 5:15 PM
 */

#ifndef _QMCALCULATOR_H
#define	_QMCALCULATOR_H

#include "qmtopology.h"

/// the idea of this class is to make QMApplications more flexible

class QMCalculator{
public:
    QMCalculator() {};
    virtual ~QMCalculator() {};

    virtual void Initialize(QMTopology *top, Property *options) {}
    virtual bool EvaluateFrame(QMTopology *top) { return true; }
    virtual void EndEvaluate(QMTopology *top) {}
protected:
};

#endif	/* _QMCALCULATOR_H */

