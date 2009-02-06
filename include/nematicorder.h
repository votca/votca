// 
// File:   nematicorder.h
// Author: ruehle
//
// Created on March 6, 2008, 4:02 PM
//

#ifndef _NEMATICORDER_H
#define	_NEMATICORDER_H

#include "topology.h"
#include "topology.h"
#include <tools/matrix.h>

class NematicOrder
{
public:
    NematicOrder() {}
    ~NematicOrder() {}
    
    void Process(Topology &top);
    
    matrix::eigensystem_t &NematicU() {return _nemat_u; }
    matrix::eigensystem_t &NematicV() {return _nemat_v; }
    matrix::eigensystem_t &NematicW() {return _nemat_w; }

private:
    matrix _mu,_mv,_mw;
    matrix::eigensystem_t _nemat_u, _nemat_v, _nemat_w;
};

#endif	/* _NEMATICORDER_H */

