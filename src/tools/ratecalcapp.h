/* 
 * File:   ratecalculator.h
 * Author: vehoff
 *
 * Created on February 2, 2010, 11:23 AM
 */

#ifndef _RATECALCULATOR_H
#define	_RATECALCULATOR_H

#include "qmapplication.h"

class RateCalculator : public QMApplication
{
public:
    RateCalculator();
    ~RateCalculator();
    
    void HelpText();
    void Initialize();
    void CheckInput();
    bool EvaluateFrame();
    void EndEvaluate();
private:
    vec _E; /// the electric field
    double _kT; /// the thermal energy

};

#endif	/* _RATECALCULATOR_H */

