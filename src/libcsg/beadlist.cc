/* 
 * File:   beadlist.cc
 * Author: ruehle
 *
 * Created on March 4, 2009, 5:42 PM
 */

#include "beadlist.h"
#include "topology.h"

int BeadList::Generate(Topology &top, const string &select)
{
    BeadContainer::iterator iter;
    
    for(iter=top.Beads().begin(); iter!=top.Beads().end();++iter) {
        if((*iter)->getType()->getName() == select) {
            push_back(*iter);
        }
    }
    return size();
}