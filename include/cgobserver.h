// 
// File:   cgobserver.h
// Author: victor
//
// Created on 4. Juni 2008, 17:08
//

#ifndef _CGOBSERVER_H
#define	_CGOBSERVER_H

#include "topology.h"

class CGObserver
{
public:
    virtual void BeginCG(Topology *top, Topology *top_atom = 0) = 0;
    virtual void EndCG() = 0;
    
    virtual void EvalConfiguration(Topology *top, Topology *top_atom = 0) = 0;
};


#endif	/* _CGOBSERVER_H */

