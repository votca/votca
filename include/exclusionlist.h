/* 
 * File:   exclusionlist.h
 * Author: ruehle
 *
 * Created on July 16, 2007, 11:10 AM
 */

#ifndef _exclusionlist_H
#define	_exclusionlist_H

#include <iostream>
#include <list>
#include <map>

using namespace std;

/// \todo fill _excl_by_bead

class ExclusionList
{
public:
    ExclusionList() {}
    ~ExclusionList() { Clear(); }
    
    void Clear(void);
    void ExcludeAll(int N);    
    void Remove(list<int> l);

    struct exclusion_t {
        int _atom;
        list<int> _exclude;
    };    
    
    //void CreateExclusions(Topology top);
    exclusion_t *GetExclusions(int bead);
    
    typedef list< exclusion_t * >::iterator iterator;
    
    iterator begin() { return _exclusions.begin(); }
    iterator end() { return _exclusions.end(); }

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

