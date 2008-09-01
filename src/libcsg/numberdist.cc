/// \addtogroup csg
///@{
// 
// File:   numberdensity.h
// Author: ruehle
//
// Created on June 9, 2008, 3:00 PM
//

#include <iostream>
#include "numberdist.h"
#include "topology.h"

NumberDist::~NumberDist(){
    Cleanup();
}

// frees allocated memory
void NumberDist::Cleanup(){
}

void NumberDist::clear() {
    for(vector<int>::iterator i = _dist.begin(); i != _dist.end(); ++i)
        *i=0;   
}

void NumberDist::Process(Configuration& conf) {
    
    double h = _cutoff / (double)_dist.size();
    for(int i = 0; i < conf.getTopology()->BeadCount(); i++){
        for(int j=i+1; j<conf.getTopology()->BeadCount(); j++){
            vec r = conf.getDist(i, j);
            double absr = abs(r);
            if(absr < _cutoff) {
		_dist[(int)(absr/h)]++;
            }
        }
    } 
}

ostream& operator<<(ostream& out, NumberDist &nd)
{
    vector<int>::iterator iter;
    double h = nd._cutoff / (double)nd._dist.size();
    double x = 0;
    for(iter=nd._dist.begin();iter!=nd._dist.end();++iter) {
        out << x << "    " << *iter << endl;
	x+=h;
    }
    return out;
}

/// @}
