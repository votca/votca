/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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
#include <fstream>
#include <boost/lexical_cast.hpp>
#include "grotopologyreader.h"


namespace votca { namespace csg {

bool GROTopologyReader::ReadTopology(string file, Topology &top)
{ 
    // cleanup topology to store new data
    ifstream fl;
    string tmp;
    top.Cleanup();
    
    fl.open(file.c_str());
    if(!fl.is_open())
        return false;

    string title;
    getline(fl, title);
    
    getline(fl, tmp);
    int natoms = atoi(tmp.c_str());
    for(;natoms>0; natoms--) {
        char c[6];
        fl.read(c, 5);
        c[5] = 0;  
    
        string resname;
        fl >> resname;
        int resnr = atoi(c);
        if(resnr > top.ResidueCount()) {
            top.CreateResidue(resname);
//            cout << " created residue " << resnr << resname<<"-\n";
        }
        string atomname;
        string x, y, z;
        fl >> atomname;
        fl >> tmp;
        fl >> x; fl >> y; fl >> z;
        top.CreateBead(1, atomname, top.GetOrCreateBeadType("no"), resnr, 1., 0.);
        getline(fl, tmp);
    }
    fl.close();
    
    cout << "WARNING: topology created from .gro file, masses and charges are wrong!\n";
    
    return true;
}

}}
