// 
// File:   boltzmanninversion.h
// Author: victor
//
// Created on 4. Juni 2008, 17:39
//

#ifndef _BONDEDSTATISTICS_H
#define	_BONDEDSTATISTICS_H

#include "cgobserver.h"
#include <tools/datacollection.h>

class BondedStatistics
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom = 0);
    void EndCG();
    
    void EvalConfiguration(Configuration *conf, Configuration *conf_atom = 0);
    
    DataCollection<double> &BondedValues() { return _bonded_values; }

protected:
    DataCollection<double> _bonded_values;
};

#endif	/* _BOLZMANNINVERSION_H */

