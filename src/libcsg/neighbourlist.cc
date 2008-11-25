// 
// File:   neighbourlist.cc
// Author: vehoff
//
// Created on March 31, 2008, 3:28 PM
//

#include <iostream>
#include "neighbourlist.h"
#include "topology.h"

NeighbourList::~NeighbourList(){
    Cleanup();
}

void NeighbourList::create_emptyNbl(int nn){
   for(int i = 0; i < nn; i++){
       entry_t* entry = new entry_t;
       _nb_list.push_back(entry);
   }
}

void NeighbourList::create_Nbs(int nb1, int nb2, vec _r){
    neighbour_t nb;
    nb._bead = (nb2);
    nb._r = _r;
    nb._dist = abs(_r);
    _nb_list[(nb1)]->_neighbours.push_back(nb);
}

void NeighbourList::Generate(Configuration& conf) {
    Topology *top = conf.getTopology();

    for(int i = 0; i < top->BeadCount(); i++){
        entry_t* entry = new entry_t;
        // creates a vector of pointers (one for each bead)
        for(int j=0; j<top->BeadCount(); j++){
            if (i!=j){ // prevents listing i as neigbour of itself
                neighbour_t nb;
                nb._r = conf.getDist(i, j);
                if((nb._dist=abs(nb._r)) < _cutoff){
                    nb._bead = j;
                    entry->_neighbours.push_back(nb);
                }
            }
        }
        // creates a list of neighbours for each molecule (and therefore each pointer)
        _nb_list.push_back(entry); // creates entire neighbour list
    } 
}

/*
// check conditions for use of gmx algorithm for triclinic pbc
void NeighbourList::CheckInput(double _cutoff, const matrix& box){
    vec a = box.getCol(0); vec b = box.getCol(1); vec c = box.getCol(2);
    if(_cutoff > 0.5*abs(c) || _cutoff > 0.5*abs(b) || _cutoff > 0.5*abs(a))
        cerr << "Cut-off too large. Violates minimum image convention. \n";
    if(a.getY() != 0 || a.getZ() != 0 || b.getZ() != 0)
        cerr << "Equation (3.1) from Gromacs Manual not fulfilled. Check your box. \n"
             << a.getY() << " " << a.getZ() << " " << b.getZ() << "\n";
    if(a.getX() <= 0 || b.getY() <= 0 || c.getZ() <= 0)
        cerr << "Equation (3.2) from Gromacs Manual not fulfilled. Check your box. \n";
    if(b.getX() > 0.5*a.getX() || b.getX() < -0.5*a.getX())
        cerr << "Equation (3.3) from Gromacs Manual not fulfilled. Check your box. \n";
    if(c.getX() > 0.5*a.getX() || c.getX() < -0.5*a.getX())
        cerr << "Equation (3.3) from Gromacs Manual not fulfilled. Check your box. \n";
    if(c.getY() > 0.5*b.getY() || c.getY() < -0.5*b.getY())
        cerr << "Equation (3.3) from Gromacs Manual not fulfilled. Check your box. \n";
}*/

// outputs the neighbor list
ostream& operator<<(ostream& out, NeighbourList &nbl)
{
   out << "There are " << nbl._nb_list.size() << " beads.\n";
   
   //vector < NeighbourList::entry_t* > 
   
   vector < NeighbourList::entry_t* >::iterator iter;
   
   for(int i=0; i<nbl._nb_list.size(); ++i){
       out << i << ":" << endl;
       NeighbourList::container::iterator iter;
       for(iter=nbl._nb_list[i]->_neighbours.begin(); iter!=nbl._nb_list[i]->_neighbours.end(); ++iter){
           out << (*iter)._bead << " r_ij: " << (*iter)._r << " d: " << (*iter)._dist << endl; 
       // (*iter) corresponds to the structure neighbour_t
       }
       out << endl;
   }
   return out;
}

// frees allocated memory
void NeighbourList::Cleanup(){
    vector<entry_t*>::iterator iter;
    for(iter=_nb_list.begin();iter!=_nb_list.end(); ++iter)
        delete (*iter);
    _nb_list.clear();
}
