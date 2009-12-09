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

