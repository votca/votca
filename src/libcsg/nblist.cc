/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include "nblist.h"
#include <iostream>

void NBList::Generate(BeadList &list1, BeadList &list2, bool do_exclusions)
{
    BeadList::iterator iter1;
    BeadList::iterator iter2;
    _do_exclusions = do_exclusions;
    
    if(list1.empty()) return;
    if(list2.empty()) return;
    
    assert(list1.getTopology() == list2.getTopology());
    Topology *top = list1.getTopology();
        
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

bool NBList::Match(Bead *bead1, Bead *bead2, const vec &r)
{
    return abs(r) < _cutoff;    
}
