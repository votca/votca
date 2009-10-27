/* 
 * File:   nblistgrid.cc
 * Author: ruehle
 *
 * Created on October 26, 2009, 4:22 PM
 */

#include <nblistgrid.h>

void NBListGrid::Generate(BeadList &list1, BeadList &list2, bool do_exclusions)
{
    BeadList::iterator iter1;
    BeadList::iterator iter2;
    _do_exclusions = do_exclusions;

    if(list1.empty()) return;
    if(list2.empty()) return;

    assert(list1.getTopology() == list2.getTopology());
    Topology *top = list1.getTopology();

    InitializeGrid(top->getBox());
    
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
            if(_do_exclusions)
                if(top->getExclusions().IsExcluded((*iter1)->getId(), (*iter2)->getId())) {
                    continue;
                }

            if(Match(*iter1, *iter2, r))
                if(!FindPair(*iter1, *iter2))
                    AddPair(new BeadPair(*iter1, *iter2, r));
        }
    }
}