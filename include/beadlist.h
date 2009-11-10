/* 
 * File:   beadlist.h
 * Author: ruehle
 *
 * Created on July 17, 2008, 5:14 PM
 */

#ifndef _BEADLIST_H
#define	_BEADLIST_H

#include <string>
#include <list>
#include "topology.h"

using namespace std;

/**
    \brief Generate lists of beads

    This class generates a list of beads based on some criteria, currently
    only the bead type.

*/

class BeadList
    : public list<Bead *>
{
public:
    BeadList() {};
    ~BeadList() {}
    
    /// \brief Select all beads of type <select>
    int Generate(Topology &top, const string &select);
    
    Topology *getTopology() {return _topology; }
    
private:
    Topology *_topology;
    
};

#endif	/* _BEADLIST_H */

