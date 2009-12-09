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

#ifndef _exclusionlist_H
#define	_exclusionlist_H

#include <iostream>
#include <list>
#include <map>

using namespace std;

/// \todo fill _excl_by_bead
/// \todo no ids but pointers, use PairList

class Topology;

class ExclusionList
{
public:
    ExclusionList() {}
    ~ExclusionList() { Clear(); }
    
    void Clear(void);
    void ExcludeAll(int N);    
    void Remove(list<int> l);
    void ExcludeList(list<int> l);
    
    struct exclusion_t {
        int _atom;
        list<int> _exclude;
    };    

    void CreateExclusions(Topology *top);
    exclusion_t *GetExclusions(int bead);    
    
    
    typedef list< exclusion_t * >::iterator iterator;
    
    iterator begin() { return _exclusions.begin(); }
    iterator end() { return _exclusions.end(); }
    
    bool IsExcluded(int bead1, int bead2);

private:
    list< exclusion_t * > _exclusions;
    map<int, exclusion_t *> _excl_by_bead;
    
    friend std::ostream &operator<<(std::ostream &out, ExclusionList& exl);
};

inline ExclusionList::exclusion_t * ExclusionList::GetExclusions(int bead)
{
   map<int, exclusion_t *>::iterator iter  = _excl_by_bead.find(bead);
   if(iter == _excl_by_bead.end()) return NULL;
   return (*iter).second;
}

std::ostream &operator<<(std::ostream &out,ExclusionList& ex);

#endif	/* _exclusionlist_H */

