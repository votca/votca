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

#include <algorithm>
#include "topology.h"
#include "exclusionlist.h"

namespace votca { namespace csg {

void ExclusionList::Clear(void)
{
    list< exclusion_t *>::iterator iter;
    
    for(iter=_exclusions.begin();iter!=_exclusions.end();++iter)
        delete *iter;    
    _exclusions.clear(); 
}

void ExclusionList::ExcludeAll(int N)
{
    Clear();
    for(int i=0; i<N-1; i++) {
        exclusion_t *e = new exclusion_t;
        e->_atom = i;
        for(int j=i+1; j<N; j++) {
            e->_exclude.push_back(j);
        }
        _exclusions.push_back(e);
    }
}

void ExclusionList::Remove(list<int> l)
{
    l.sort();
    list<int>::iterator i, j;
    list<exclusion_t*>::iterator ex;
    
    for(i=l.begin(); i!=l.end(); i++) {
        for(ex=_exclusions.begin(); ex!=_exclusions.end(); ++ex)
            if((*ex)->_atom == *i) break;
        if(ex==_exclusions.end()) continue;
        j = i;        
        for(++j; j!=l.end(); j++)
            (*ex)->_exclude.remove(*j);
        
        if((*ex)->_exclude.empty()) {
        //    delete *ex;
            (*ex)=NULL;
            _exclusions.erase(ex);
        }
    }
    _exclusions.remove(NULL);
}

void ExclusionList::ExcludeList( list<int> l ) {
    l.sort();
    list<int>::iterator i, j, k;
    list<exclusion_t*>::iterator ex;
    
    for ( i = l.begin(); i != l.end(); ++i ) {
        for (ex = _exclusions.begin(); ex != _exclusions.end(); ++ex)
            if ( (*ex)->_atom == (*i) ) break;
        if (ex==_exclusions.end()) { // so far there are no exclusions for i
            exclusion_t *e = new exclusion_t;
            e->_atom = (*i);
            
            j = i;
            for (++j; j != l.end(); ++j)
                e->_exclude.push_back( (*j) );
            _exclusions.push_back(e);
              
        }
        else {
            // there are some exclusions for i already. Add new exclusions if they are 
            // not there yet!
            j = i;
            for (++j; j != l.end(); ++j) {
                
                for ( k = (*ex)->_exclude.begin(); k != (*ex)->_exclude.end(); ++k ) 
                    if ( (*j) == (*k) ) break;
                if ( k == (*ex)->_exclude.end() ) (*ex)->_exclude.push_back( (*j) );
            }
            
        }
    }
    
}
    
void ExclusionList::CreateExclusions(Topology *top) {
    
    InteractionContainer &ic = top->BondedInteractions();
    InteractionContainer::iterator ia;
    list<int> l;
    
    for (ia = ic.begin(); ia != ic.end(); ++ia) {
        int beads_in_int = (*ia)->BeadCount();
        l.clear();
        
        for (int ibead = 0; ibead < beads_in_int; ibead ++) {
            int ii = (*ia)->getBeadId(ibead);
            l.push_back(ii);
        }
        ExcludeList(l);  
    }
    // Create map
    list< exclusion_t * >::iterator it;
    for ( it = _exclusions.begin(); it != _exclusions.end(); ++it ) {
        int iatom = (*it)->_atom;
        _excl_by_bead[ iatom ] = (*it);
    }
    
}

bool ExclusionList::IsExcluded(int bead1, int bead2) {
    exclusion_t *excl;
    if (bead2 < bead1) swap(bead1, bead2);
    if ((excl = GetExclusions(bead1))) {
        if(find(excl->_exclude.begin(), excl->_exclude.end(), bead2) 
                != excl->_exclude.end()) return true;        
    }
    return false;
}

std::ostream &operator<<(std::ostream &out, ExclusionList& exl)
{
    list<ExclusionList::exclusion_t*>::iterator ex;
    
    for(ex=exl._exclusions.begin();ex!=exl._exclusions.end();++ex) {
        list<int>::iterator i;
        out << (int)((*ex)->_atom) + 1;
        for(i=(*ex)->_exclude.begin(); i!=(*ex)->_exclude.end(); ++i) {
            out << " " << ((*i)+1);
        }
        out << endl;
    }
    return out;
}

void ExclusionList::InsertExclusion(int index, list<int> l) {

    list<int>::iterator i;

    //cout << "atom " << index << endl;
    for(i=l.begin(); i!=l.end(); ++i) {
        int bead1 = index;
        int bead2 = *i;
        if (bead2 < bead1) swap(bead1, bead2);
        exclusion_t *e;
        if((e = GetExclusions(bead1)) == NULL) {
            e = new exclusion_t;
            e->_atom = bead1;
            _exclusions.push_back(e);
        }
        e->_exclude.push_back(bead2);
        //cout << (*i) << " " << endl;
    }
    //cout << endl;
}

}}
