/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_CSG_EXCLUSIONLIST_H
#define	_VOTCA_CSG_EXCLUSIONLIST_H

#include <iostream>
#include <list>
#include <map>
#include "bead.h"

namespace votca { namespace csg {
using namespace votca::tools;

/// \todo fill _excl_by_bead
/// \todo no ids but pointers, use PairList

class Topology;
class Bead;

class ExclusionList
{
public:
    ExclusionList() {}
    ~ExclusionList() { Clear(); }
    
    void Clear(void);

    template<typename iteratable>
    void Remove(iteratable &l);

    template<typename iteratable>
	void ExcludeList(iteratable &l);
    
    struct exclusion_t {
        Bead *_atom;
        std::list<Bead *> _exclude;
    };    

    void CreateExclusions(Topology *top);
    exclusion_t *GetExclusions(Bead *bead);
    
    typedef std::list< exclusion_t * >::iterator iterator;
    
    iterator begin() { return _exclusions.begin(); }
    iterator end() { return _exclusions.end(); }
    
    bool IsExcluded(Bead *bead1, Bead *bead2);

    template<typename iteratable>
    void InsertExclusion(Bead *bead, iteratable &excluded);

    void InsertExclusion(Bead *bead1, Bead *bead2);

    void RemoveExclusion(Bead *bead1, Bead *bead2);
private:
    std::list< exclusion_t * > _exclusions;
    std::map<Bead *, exclusion_t *> _excl_by_bead;
    
    friend std::ostream &operator<<(std::ostream &out, ExclusionList& exl);
};

inline ExclusionList::exclusion_t * ExclusionList::GetExclusions(Bead *bead)
{
  std::map<Bead *, exclusion_t *>::iterator iter  = _excl_by_bead.find(bead);
   if(iter == _excl_by_bead.end()) return NULL;
   return (*iter).second;
}

template<typename iteratable>
inline void ExclusionList::Remove(iteratable &l)
{
    typename iteratable::iterator i, j;

    for ( i = l.begin(); i != l.end(); ++i ) {
    	for ( j = i; j != l.end(); ++j ) {
    		RemoveExclusion(*i, *j);
    	}
    }
}

template<typename iteratable>
inline void ExclusionList::ExcludeList( iteratable &l ) {
    typename iteratable::iterator i, j;

    for ( i = l.begin(); i != l.end(); ++i ) {
    	for ( j = i; j != l.end(); ++j ) {
    		InsertExclusion(*i, *j);
    	}
    }
}

template<typename iteratable>
inline void ExclusionList::InsertExclusion(Bead *bead1_, iteratable &l)
{
	for(typename iteratable::iterator i=l.begin(); i!=l.end(); ++i) {
		Bead *bead1 = bead1_;
		;Bead *bead2 = *i;
		if (bead2->getId() < bead1->getId()) std::swap(bead1, bead2);
		if(bead1==bead2) continue;
		if(IsExcluded(bead1, bead2)) continue;
		exclusion_t *e;
		if((e = GetExclusions(bead1)) == NULL) {
			e = new exclusion_t;
			e->_atom = bead1;
			_exclusions.push_back(e);
			_excl_by_bead[ bead1 ] = e;
		}
		e->_exclude.push_back(bead2);
	}
}

//template<>
inline void ExclusionList::InsertExclusion(Bead *bead1, Bead *bead2) {
    if (bead2->getId() < bead1->getId()) std::swap(bead1, bead2);
	if(bead1==bead2) return;
	if(IsExcluded(bead1, bead2)) return;

	exclusion_t *e;
	if((e = GetExclusions(bead1)) == NULL) {
		e = new exclusion_t;
		e->_atom = bead1;
		_exclusions.push_back(e);
		_excl_by_bead[ bead1 ] = e;
	}
	e->_exclude.push_back(bead2);
}

inline void ExclusionList::RemoveExclusion(Bead *bead1, Bead *bead2) {
    if (bead2->getId() < bead1->getId()) std::swap(bead1, bead2);
    if(bead1==bead2) return;
    if(!IsExcluded(bead1, bead2)) return;
    std::list<exclusion_t*>::iterator ex;
    for(ex=_exclusions.begin(); ex!=_exclusions.end(); ++ex)
        if((*ex)->_atom == bead1) break;
    if(ex==_exclusions.end()) return;
    (*ex)->_exclude.remove(bead2);
    if((*ex)->_exclude.empty()) {
        (*ex)=NULL;
        _exclusions.erase(ex);
    }
    _exclusions.remove(NULL);
}

std::ostream &operator<<(std::ostream &out,ExclusionList& ex);

}}

#endif	/* _VOTCA_CSG_EXCLUSIONLIST_H */



