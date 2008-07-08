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

void NeighbourList::Generate(Configuration& conf) {
    CheckInput(_cutoff, conf.getBox());
    for(int i = 0; i < conf.getTopology()->BeadCount(); i++){
        entry_t* entry = new entry_t;
        // creates a vector of pointers (one for each bead)
        for(int j=0; j<conf.getTopology()->BeadCount(); j++){
            if (i!=j){ // prevents listing i as neigbour of itself
                neighbour_t nb;
                nb._r = CalcDist(conf.getPos(i), conf.getPos(j), conf.getBox());
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

// calculates closest distance between mol A and mol B or pbc images thereof using
// minimum image convention and a box satisfying conditions checked by CheckInput
vec NeighbourList::CalcDist(const vec& r_i, const vec& r_j, const matrix& box){
    vec r_tp, r_dp, r_sp, r_ij;
    vec a = box.getCol(0); vec b = box.getCol(1); vec c = box.getCol(2);
    r_tp = r_j - r_i;
    r_dp = r_tp - c*round(r_tp.getZ()/c.getZ());  
    r_sp = r_dp - b*round(r_dp.getY()/b.getY());
    r_ij = r_sp - a*round(r_sp.getX()/a.getX());
    #ifdef DEBUG
    cout << "r_ij: " << r_ij << " r_tp: " << r_tp << " d: " << abs(r_ij) << "\n";
    #endif
    return r_ij;
}

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
}

void NeighbourList::getPbc(Configuration& conf){
     matrix box = conf.getBox();
     vec a = box.getCol(0); vec b = box.getCol(1); vec c = box.getCol(2);
     cout << "The box vectors are [nm]:" << endl;
     cout << "a= " << a << endl << "b= " << b << endl << "c= " << c << endl;
     cout << "Cut-off: " << getCutoff() << " [nm]" << endl;
}

// outputs the neighbor list
ostream& operator<<(ostream& out, NeighbourList &nbl)
{
   out << "There are " << nbl._nb_list.size() << " beads.\n";
   
   // vector<NeighbourList::entry_t*>::iterator iter;
   for(int i=0; i<nbl._nb_list.size(); ++i){
       out << i << ":" << endl; 
       NeighbourList::container::iterator iter;
       for(iter=nbl._nb_list[i]->_neighbours.begin(); iter!=nbl._nb_list[i]->_neighbours.end(); iter++){
           out << (*iter)._bead << " r_ij: " << (*iter)._r << " d: " << (*iter)._dist << endl; 
       // (*iter) corresponds to the structure neighbour_t
       }
       out << "\n";
   }
   /* int lastmol = nbl._nb_list.size()-1;
   out << "Vectors from last molecule (" << lastmol << ") to neighbors: \n";
   NeighbourList::container::iterator iter;
       for(iter=nbl._nb_list[lastmol]->_neighbours.begin(); iter!=nbl._nb_list[lastmol]->_neighbours.end(); iter++){
           out << (*iter)._bead << ": " << (*iter)._r << "\n"; 
       }
    */
   return out;
}

// frees allocated memory
void NeighbourList::Cleanup(){
    vector<entry_t*>::iterator iter;
    for(iter=_nb_list.begin();iter!=_nb_list.end(); ++iter)
        delete (*iter);
    _nb_list.clear();
}
