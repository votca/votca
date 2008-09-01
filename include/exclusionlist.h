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

using namespace std;

class ExclusionList
{
public:
    ExclusionList() {}
    ~ExclusionList() { Clear(); }
    
    void Clear(void);
    void ExcludeAll(int N);    
    void Remove(list<int> l);

private:
    struct exclusion_t {
        int _atom;
        list<int> _exclude;
    };
    list< exclusion_t * > _exclusions;
    
    friend std::ostream &operator<<(std::ostream &out, ExclusionList& exl);
};

std::ostream &operator<<(std::ostream &out,ExclusionList& ex);

#endif	/* _exclusionlist_H */

