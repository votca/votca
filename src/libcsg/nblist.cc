/* 
 * File:   nblist.cc
 * Author: ruehle
 *
 * Created on March 6, 2009, 2:24 PM
 */

#include "nblist.h"
#include <iostream>

void NBList::Generate(BeadList &list1, BeadList &list2, ExclusionList *ExcList=0)
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
            if(Match(*iter1, *iter2, r, ExcList))
                if(!FindPair(*iter1, *iter2))
                    AddPair(new BeadPair(*iter1, *iter2, r));
        } 
    }
}

bool NBList::Match(Bead *bead1, Bead *bead2, const vec &r, ExclusionList *ExcList=0)
{

    int bead1_id = bead1->getId();
    int bead2_id = bead2->getId();  
    
    if (bead1_id > bead2_id) {
        int tmp = bead1_id;
        bead1_id = bead2_id;
        bead2_id = tmp;
    }
    ExclusionList::exclusion_t *excl;
    if (ExcList != 0) {
        if(excl = ExcList->GetExclusions(bead1_id)) {
            list<int> exclusions_bead1 = excl->_exclude;   
            list<int>::iterator iter;
            cout << "excl_size=" << exclusions_bead1.size() << endl;
    
            for (iter = exclusions_bead1.begin(); iter != exclusions_bead1.end(); ++iter)
                if ((*iter)==bead2_id) return false;
        }
    }
    return abs(r) < _cutoff;    
}
