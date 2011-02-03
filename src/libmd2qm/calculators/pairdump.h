/* 
 * File:   calc_integrals.h
 * Author: vehoff
 *
 * Created on April 8, 2010, 10:33 AM
 */

#ifndef _PAIR_DUMP_H
#define	_PAIR_DUMP_H

#include "qmpair.h"
#include "paircalculator.h"
#include <sys/stat.h>

class PairDump : public PairCalculator
{
public:
    PairDump() {};
    ~PairDump() {};

    void EvaluatePair(QMTopology *top, QMPair *pair);
};

inline void PairDump::EvaluatePair(QMTopology *top, QMPair *pair){
    CrgUnit *crg1 = pair->Crg1();
    CrgUnit *crg2 = pair->Crg2();

    string framedir=string("frame")+lexical_cast<string>(top->getStep()) +string("/") ;
    mkdir(framedir.c_str(),0755);
    top->GetJCalc().WriteProJ(*crg1, *crg2,framedir);
        
}

#endif	
