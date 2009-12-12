/* 
 * File:   qmtopology.h
 * Author: james
 *
 * Created on 12 December 2009, 13:14
 */

#ifndef _ATOMQMTOPOLOGY_H
#define	_ATOMQMTOPOLOGY_H

#include <votca/csg/topology.h>
#include <votca/csg/nblist.h>

/**
    \brief contains the atomistic representation of the crg units

 * It contains Beads describing atoms. It should contain the overloaded CreateBeads function which given
 * a crgunit and a chrage unit type generates all the relevant atomistic beads
*/

class AtomQMTopology : Topology{
public:
    QMTopology();
    ~QMTopology();
    int CreateBeads();
    int Init();
private:

};

#endif	/* _QMTOPOLOGY_H */

