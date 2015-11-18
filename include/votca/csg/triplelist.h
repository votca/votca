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

#ifndef _TRIPLELIST_H
#define	_TRIPLELIST_H

#include <list>
#include <map>

namespace votca { namespace csg {
using namespace std;

template<typename element_type, typename triple_type>
class TripleList {
public:
    TripleList() {}
    virtual ~TripleList() { Cleanup(); }
           
    void AddTriple(triple_type *t);
    
    typedef typename std::list<triple_type *>::iterator iterator;
    
    iterator begin() { return _triples.begin(); }
    iterator end() { return _triples.end(); }
    typename list<triple_type*>::size_type size() { return _triples.size(); }    
    triple_type *front() { return _triples.front(); }
    triple_type *back() { return _triples.back(); }    
    bool empty() { return _triples.empty(); }
    
    void Cleanup();
    
    triple_type *FindTriple(element_type e1, element_type e2, element_type e3);

    typedef element_type element_t;
    typedef triple_type triple_t;

protected:
    list<triple_type *> _triples;
      
    map< element_type , map<element_type, map<element_type, triple_type *> > > _triple_map;
    
};

template<typename element_type, typename triple_type>
inline void TripleList<element_type, triple_type>::AddTriple(triple_type *t)
{
    /// \todo be careful, same triple object is used, some values might change (e.g. sign of distance vectors)
    //experimental: So far only mapping '123' and '321' to the same triple
    _triple_map[ (*t)[0] ][ (*t)[1] ][ (*t)[2] ] = t;
    _triple_map[ (*t)[2] ][ (*t)[1] ][ (*t)[0] ] = t;
     /// \todo check if unique    
    _triples.push_back(t);    
}


template<typename element_type, typename triple_type>
inline void TripleList<element_type, triple_type>::Cleanup()
{
    for(iterator iter = _triples.begin(); iter!=_triples.end(); ++iter)
        delete *iter;
    _triples.clear();
    _triple_map.clear();
}

template<typename element_type, typename triple_type>
inline triple_type *TripleList<element_type, triple_type>::FindTriple(element_type e1, element_type e2, element_type e3)
{
    typename std::map< element_type , map<element_type, map<element_type, triple_type *> > > ::iterator iter1;    
    iter1 = _triple_map.find(e1);
    if(iter1==_triple_map.end()) return NULL;
     
    typename std::map<element_type, map<element_type, triple_type *> >::iterator iter2;    
    iter2 = iter1->second.find(e2);
    if(iter2 == iter1->second.end()) return NULL;
    
    typename std::map<element_type, triple_type *>::iterator iter3;    
    iter3 = iter2->second.find(e3);
    if(iter3 == iter2->second.end()) return NULL;
    
    return iter3->second;
}

}}

#endif	/* _TRIPLELIST_H */

