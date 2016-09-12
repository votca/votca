/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/csg/nblist_3body.h>
#include <iostream>

namespace votca { namespace csg {

NBList_3Body::NBList_3Body()
  : _do_exclusions(false), _match_function(0)
{
    setTripleType<BeadTriple>();
    SetMatchFunction(NBList_3Body::match_always);
}

NBList_3Body::~NBList_3Body()
{
  // TODO: NBList_3Body destructor
}

void NBList_3Body::Generate(BeadList &list1, BeadList &list2, BeadList &list3, bool do_exclusions)
{
    BeadList::iterator iter1;
    BeadList::iterator iter2;
    BeadList::iterator iter3;
    _do_exclusions = do_exclusions;
    
    if(list1.empty()) return;
    if(list2.empty()) return;
    if(list3.empty()) return;
    
    assert(list1.getTopology() == list2.getTopology());
    assert(list1.getTopology() == list3.getTopology());
    assert(list2.getTopology() == list3.getTopology());
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
            if(&list2 == &list3) {        
                iter3=iter2;
                ++iter3;
            }
            if(&list1 == &list3) {        
                iter3=iter1;
                ++iter3;
            }
            if( (&list2 != &list3) && (&list1 != &list3) ) {
                iter3 = list3.begin();                
            }            
            
            if(*iter1 == *iter3) continue;
            if(*iter2 == *iter3) continue;
            
            for(; iter3 != list3.end(); ++iter3) {
                vec u = (*iter1)->getPos();
                vec v = (*iter2)->getPos();
                vec z = (*iter3)->getPos();
            
                vec r12 = top->BCShortestConnection(u, v);
                vec r13 = top->BCShortestConnection(u, z);                
                vec r23 = top->BCShortestConnection(v, z);                
                double d12 = abs(r12);
                double d13 = abs(r13);                
                double d23 = abs(r23);
                //to do: at the moment use only one cutoff value
                if( (d12 < _cutoff) && (d13 < _cutoff) && (d23 < _cutoff) ){
                /// experimental: at the moment exclude interaction as soon as one of the three pairs (1,2) (1,3) (2,3) is excluded!
                if(_do_exclusions)
                    if( (top->getExclusions().IsExcluded(*iter1, *iter2)) || (top->getExclusions().IsExcluded(*iter1, *iter3)) || (top->getExclusions().IsExcluded(*iter2, *iter3)) ) {
                        continue;
                    }
                    if((*_match_function)(*iter1, *iter2, *iter3, r12, r13, r23, d12, d13, d23))
                        if(!FindTriple(*iter1, *iter2, *iter3))
                            AddTriple( _triple_creator(*iter1, *iter2, *iter3, r12, r13, r23));
                }
                
            }
        } 
    }    

}

}}
