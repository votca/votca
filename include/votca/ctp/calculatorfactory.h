/* 
 * File:   CalculatorFactory.h
 * Author: ruehle
 *
 * Created on June 30, 2010, 4:32 PM
 */

#ifndef _CALCULATORFACTORY_H
#define	_CALCULATORFACTORY_H

#include <votca/tools/objectfactory.h>
#include "qmcalculator.h"

using namespace std;

class CalculatorFactory
: public ObjectFactory<std::string, QMCalculator>
{
private:
    CalculatorFactory() {}
public:
    
    static void RegisterAll(void);

    friend CalculatorFactory &Calculators();
};

inline CalculatorFactory &Calculators()
{
    static CalculatorFactory _instance;
    return _instance;
}

#endif	/* _CALCULATORFACTORY_H */

