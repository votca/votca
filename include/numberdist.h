// 
// File:   numberdist.h
// Author: ruehle
//
// Created on June 9, 2008, 3:00 PM
//

#ifndef _NUMBERDIST_H
#define	_NUMBERDIST_H

#include <list>
#include <vector>
#include "topology.h"

using namespace std;

/**
    Calculates the number distrubution N(r) (or S(r)) 
*/
class NumberDist
{
public:
     NumberDist() {}; // constructor
     ~NumberDist(); // destructor
     
     void Cleanup(); // clear occupied memory
     
     // define cut-off radius
     void setCutoff(double cutoff) { _cutoff = cutoff; }
     // show used cut-off radius 
     double getCutoff() { return _cutoff; } 

     void setN(int N) { _dist.resize(N); }
     int getN() { return _dist.size(); }
 
     // creates entire neighbour list
     // void GenerateIntraMolecular(Configuration &conf);
     // void GenerateInterMolecular(Configuration &conf);
     void Process(Topology &top); 

     void clear();
     
     vector<int> &getDist() { return _dist; }
private:
    double _cutoff;
    
    vector<int> _dist;
        
    friend ostream& operator<<(ostream& out, NumberDist &nd);
};

// output numberdist
ostream& operator<<(ostream& out, NumberDist &nd);

#endif	/* _NUMBERDENSITY_H */

