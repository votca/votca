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


class CrgTopology:Topology{
public:
    CrgTopology();
    ~CrgTopology();

private:

    /// the list of neighbours
    NBList * _neighs;
};

#endif	/* _CRGTOPOLOGY_H */

