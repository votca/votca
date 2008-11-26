// 
// File:   neighbourlist.h
// Author: ruehle
//
// Created on March 31, 2008, 2:00 PM
//

#ifndef _NEIGHBOURLIST_H
#define	_NEIGHBOURLIST_H

#include <list>
#include <vector>
#include "configuration.h"

using namespace std;

class NeighbourList
{
public:
     NeighbourList(){}; /// constructor
     ~NeighbourList(); /// destructor
     
     void Cleanup(); /// clear occupied memory
     
     void create_emptyNbl(int nn); /// creates empty nbl for nn entries
     
     void setCutoff(double cutoff) { _cutoff = cutoff; } /// define cut-off radius     
     double getCutoff() { return _cutoff; } /// show used cut-off radius
               
     void Generate(Configuration &conf); // creates entire neighbour list
     // void GenerateIntraMolecular(Configuration &conf);
     // void GenerateInterMolecular(Configuration &conf);
     
     
     void create_Nbs(int nb1, int nb2, vec _r);
     
     // this is a neighbour
     struct neighbour_t {    
         int _bead;     // the id of the neighbour
         // maybe also store pointer to neighbour in list?
         // entry_t *_nb;
         vec _r;        // the vector from current to neighbour
         double _dist;  // the distance of the neighbours, do we need that, double information...
     };
    
     // not sure weather it's better to store as list or vector
     // this typedef makes it flexible, use of iterators makes it easy to change
     // in between these two
     // think about: is it better to store pointers, but then have to cleanup
     typedef list<neighbour_t> container; 

     // find a better name for that struct
     // this contains all the neighbours of 1 bead
     struct entry_t {
         container _neighbours;
         // do we need more here?
     };

     // allows access to vector<entry_t*>
     vector<entry_t*> &NbList() { return _nb_list; }

     /// \todo optimize this, not very clean implementation
     struct pair_t {
         int bead1, bead2;
         vec r;
         double dist;

         pair_t(int b1, int b2, vec rr, double d)
         : bead1(b1), bead2(b2), r(rr), dist(d) {}
     };
     
    struct pair_iterator {
        pair_iterator() : _nb_list(0) {}
        
        pair_t operator*() const
            { return pair_t(_curbead, (*_curneigh)._bead, (*_curneigh)._r, (*_curneigh)._dist); }
        
        pair_iterator &operator++();                      
        
        bool operator==(const pair_iterator &);                      
        bool operator!=(const pair_iterator &);                      
         
    private:
        pair_iterator(vector<entry_t*> &nb_list, bool bEndIterator=false);
        
        vector<NeighbourList::entry_t *> _nb_list;
        int _curbead;
        NeighbourList::container::iterator _curneigh;
        
        friend class NeighbourList;
    };
    
    pair_iterator begin_pairs() { return pair_iterator(_nb_list); }
    pair_iterator end_pairs() { return pair_iterator(_nb_list, true); }

private:
    // this stores the neighbours of all beads
    // one entry for each bead
    vector<entry_t *> _nb_list;
    double _cutoff;
    
    friend ostream& operator<<(ostream& out, NeighbourList &nb);
};

// output neighbourlist, look at data collection & how it works
ostream& operator<<(ostream& out, NeighbourList &nb);


inline NeighbourList::pair_iterator &NeighbourList::pair_iterator::operator++()
{
    _curneigh++;
    // did we read the last neighbour of current bead?
    if(_curneigh==_nb_list[_curbead]->_neighbours.end()) {
        // look for next valid entry, as long as there are still beads
        while(_curbead < _nb_list.size()) {
            // go through the neighbours of the bead
            for(_curneigh = _nb_list[_curbead]->_neighbours.begin();
                _curneigh!=_nb_list[_curbead]->_neighbours.end();
                ++_curneigh) {   
                // break if we found neighbour with higher id
                if((*_curneigh)._bead > _curbead)
                    return *this;
            }
            // go to the next bead
            ++_curbead;                
        }
        // no valid bead found, point to end
        _curbead = _nb_list.size();
        _curneigh = _nb_list[0]->_neighbours.begin();
        return *this;
    }    
    return *this;    
}
        
inline bool NeighbourList::pair_iterator::operator==(const NeighbourList::pair_iterator &i)
{
    return ((_curbead == i._curbead) && (_curneigh == i._curneigh));
}

bool NeighbourList::pair_iterator::operator!=(const NeighbourList::pair_iterator &i)
{
    return !((*this)==i);
}

NeighbourList::pair_iterator::pair_iterator(vector<NeighbourList::entry_t*> &nb_list, bool bEndIterator)
        : _nb_list(nb_list)
{
    if(!bEndIterator) {
        _curbead=0;
        _curneigh=_nb_list[_curbead]->_neighbours.begin();
        return;
    }
    _curbead=_nb_list.size();
    _curneigh=_nb_list[0]->_neighbours.begin();    
}

#endif	/* _NEIGHBOURLIST_H */

