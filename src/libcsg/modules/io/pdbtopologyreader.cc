/// \addtogroup csg
///@{
// 
// File:   pdbtopologyreader.cc
// Author: victor
//
// Created on 27. Dezember 2007, 21:19
//

#include <iostream>
#include "pdbtopologyreader.h"

bool PDBTopologyReader::ReadTopology(string file, Topology &top)
{ 
    // cleanup topology to store new data
    /*top.Cleanup();
    
    // use gmx functions to read topology file
    // read the residue names
    for(int i=0; i < gtp.atoms.nres; i++) {
        top.CreateResidue(*(gtp.atoms.resname[i]));
        //cout << "residue " << i << ": " << *(gtp.atoms.resname[i]) << endl;
    }
    
    // read the atoms
    for(int i=0; i < gtp.atoms.nr; i++) {
        top.CreateBead(1, *(gtp.atoms.atomname[i]), a->resnr, a->m, a->q);  
        //cout << *(gtp.atoms.atomname[i]) << " residue: " << a->resnr << endl;
    }
    
    BeadInfo *bead;
    for(int i=0; i<gtp.blocks[gmx::ebMOLS].nr; i++) {
        int i1 = gtp.blocks[gmx::ebMOLS].index[i];
        int i2 = gtp.blocks[gmx::ebMOLS].index[i+1];
        MoleculeInfo *mi = top.CreateMolecule("PPY10");    
        int res0 = top.getBead(i1)->getResnr();
        for(int i=i1; i<i2; ++i) {
            bead = top.getBead(i);
            stringstream n("");
            
            n << bead->getResnr()-res0 + 1 << ":" <<  top.getResidue(bead->getResnr())->getName() << ":" << bead->getName();
            //cout << n.str() << endl;
            mi->AddBead(i, n.str());
        }         
    }
    
    return true;/**/
}
/// @}
