/* 
 * File:   crgtopology.h
 * Author: james
 *
 * Created on 12 December 2009, 13:50
 */

#ifndef _CRGTOPOLOGY_H
#define	_CRGTOPOLOGY_H

/**
    \brief this class knows which beads make which crg unit

 
*/

#include <votca/csg/topology.h>
#include <votca/csg/nblist.h>
#include <moo/crgunittype.h>
#include <moo/crgunit.h>

#include "qmbead.h"
/**
    \brief contains the beads describing the c.o.m. of each cahrge unit
 * given a CG toplogy, it should be apply the mapping for the cg beads -> qm beads and
 * it should update the position of the crg unit. Crg units should be associated to
 * these qm beads and not to any other.
*/

class CrgTopology:Topology{
public:
    CrgTopology();
    ~CrgTopology();

private:

    /// the list of neighbours
    NBList * _neighs;
    vector <QMBead *>;
};

#endif	/* _CRGTOPOLOGY_H */

