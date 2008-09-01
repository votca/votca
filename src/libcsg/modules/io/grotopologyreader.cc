/// \addtogroup csg
///@{
// 
// File:   grotopologyreader.cc
// Author: victor
//
// Created on 27. Dezember 2007, 21:47
//

// 
// File:   pdbtopologyreader.cc
// Author: victor
//
// Created on 27. Dezember 2007, 21:19
//

#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include "grotopologyreader.h"

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
        Residue *res = top.getResidue(resnr-1);
        string atomname;
        string x, y, z;
        fl >> atomname;
        fl >> tmp;
        fl >> x; fl >> y; fl >> z;
        top.CreateBead(1, atomname, resnr, 1, 0);
        getline(fl, tmp);
    }
    fl.close();
    
    cout << "WARNING: topology created from .gro file, masses and charges are wrong!\n";
    
    return true;
}
/// @}
