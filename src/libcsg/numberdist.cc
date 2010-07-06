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

#include <iostream>
#include "numberdist.h"
#include "topology.h"

namespace votca { namespace csg {

NumberDist::~NumberDist(){
    Cleanup();
}

// frees allocated memory
void NumberDist::Cleanup(){
}

void NumberDist::clear() {
    for(vector<int>::iterator i = _dist.begin(); i != _dist.end(); ++i)
        *i=0;   
}

void NumberDist::Process(Topology& top) {
    
    double h = _cutoff / (double)_dist.size();
    for(int i = 0; i < top.BeadCount(); i++){
        for(int j=i+1; j< top.BeadCount(); j++){
            vec r = top.getDist(i, j);
            double absr = abs(r);
            if(absr < _cutoff) {
		_dist[(int)(absr/h)]++;
            }
        }
    } 
}

ostream& operator<<(ostream& out, NumberDist &nd)
{
    vector<int>::iterator iter;
    double h = nd._cutoff / (double)nd._dist.size();
    double x = 0;
    for(iter=nd._dist.begin();iter!=nd._dist.end();++iter) {
        out << x << "    " << *iter << endl;
	x+=h;
    }
    return out;
}

}}
