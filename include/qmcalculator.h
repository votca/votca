/* 
 * File:   qmcalculator.h
 * Author: vehoff
 *
 * Created on March 31, 2010, 5:15 PM
 */

#ifndef _QMCALCULATOR_H
#define	_QMCALCULATOR_H

#include <boost/program_options.hpp>
#include "qmtopology.h"
#include <votca/tools/property.h>

/// the idea of this class is to make QMApplications more flexible

class QMCalculator{
public:
    QMCalculator();
    ~QMCalculator();

    virtual void Initialize();
    virtual void EvaluateFrame();
    virtual void EndCalc();
};

#endif	/* _QMCALCULATOR_H */

