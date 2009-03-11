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

class BeadList
    : public list<Bead *>
{
public:
    BeadList() {}
    ~BeadList() {}
    
    int Generate(Topology &top, const string &select);
    
private:
};

#endif	/* _BEADLIST_H */

