/* 
 * File:   nblist.cc
 * Author: ruehle
 *
 * Created on March 6, 2009, 2:24 PM
 */

#include "nblist.h"

void NBList::Generate(BeadList &list1, BeadList &list2)
{
    BeadList::iterator iter1;
    BeadList::iterator iter2;
    
    if(list1.empty()) return;
    if(list2.empty()) return;
    
    assert(list1.front()->getParent() == list2.front()->getParent());
    
    Topology *top = list1.front()->getParent();
    
    for(iter1 = list1.begin(); iter1 != list1.end(); ++iter1) {
        if(&list1 == &list2) {        
            iter2=iter1;
            ++iter2;
        }
        else
            iter2 = list2.begin();
            
        if(*iter1 == *iter2) continue;
    
        for(; iter2 != list2.end(); ++iter2) {
            vec u = (*iter1)->getPos();
            vec v = (*iter2)->getPos();
            
            vec r = top->BCShortestConnection(u, v);
            if(Match(*iter1, *iter2, r))
                if(!FindPair(*iter1, *iter2))
                    AddPair(new BeadPair(*iter1, *iter2, r));
        } 
    }
}

bool NBList::Match(Bead *bead1, Bead *bead2, const vec &r)
{
    return abs(r) < _cutoff;
}
