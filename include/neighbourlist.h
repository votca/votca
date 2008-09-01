/// \addtogroup csg
///@{
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

// think about: do we need s.th. to treat molecules,
// e.g. only intra molecule, inter-molecule, ..?

class NeighbourList
{
public:
     NeighbourList(){}; // constructor
     ~NeighbourList(); // destructor
     
     void Cleanup(); // clear occupied memory
     
     void setCutoff(double cutoff) { _cutoff = cutoff; } // define cut-off radius
     double getCutoff() { return _cutoff; } // show used cut-off radius
     void CheckInput(double _cutoff, const matrix& box); // check box definition & cutoff
     
     vec CalcDist(const vec& r_i, const vec& r_j, const matrix &box);
     // calculates distance to surrounding neighbors including closest images
     // due to periodic boundary conditions which may be triclinic
     
     void Generate(Configuration &conf); // creates entire neighbour list
     // void GenerateIntraMolecular(Configuration &conf);
     // void GenerateInterMolecular(Configuration &conf);
     
     void getPbc(Configuration& conf); // prints box vectors and cut-off
    
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
    
private:
    // this stores the neighbours of all beads
    // one entry for each bead
    vector<entry_t *> _nb_list;
    double _cutoff;
    
    friend ostream& operator<<(ostream& out, NeighbourList &nb);
};

// output neighbourlist, look at data collection & how it works
ostream& operator<<(ostream& out, NeighbourList &nb);

#endif	/* _NEIGHBOURLIST_H */

/// @}
