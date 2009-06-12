// 
// File:   gmxtopologyreader.cc
// Author: ruehle
//
// Created on April 5, 2007, 12:35 PM
//

#include <iostream>
#include "gmxtopologyreader.h"

namespace gmx {
   extern "C" {
        #include <statutil.h>
        #include <typedefs.h>
        #include <smalloc.h>
        #include <vec.h>
        #include <copyrite.h>
        #include <statutil.h>
        #include <tpxio.h>
    }
    // this one is needed because of bool is defined in one of the headers included by gmx
    #undef bool
}

bool GMXTopologyReader::ReadTopology(string file, Topology &top)
{ 
    gmx::t_topology gtp;
    char       title[STRLEN];
    gmx::rvec       *xtop;
    gmx::matrix     box;
    int ePBC;
    // cleanup topology to store new data
    top.Cleanup();
    
    // use gmx functions to read topology file
    if(gmx::read_tps_conf((char *)file.c_str(),title, &gtp,&ePBC, &xtop,NULL,box,TRUE) == TRUE)
    gmx::sfree(xtop);
    //cout << "number of blocks: " << gtp.blocks[gmx::ebMOLS].nra << endl;
    //cout << "number of blocks: " << gtp.blocks[gmx::ebMOLS].a[0] << endl;
    // read the residue names
    for(int i=0; i < gtp.atoms.nres; i++) {
#ifdef GMX4CVS
        top.CreateResidue(*(gtp.atoms.resinfo[i].name));
#else
        top.CreateResidue(*(gtp.atoms.resname[i]));            
#endif
        //cout << "residue " << i << ": " << *(gtp.atoms.resname[i]) << endl;
    }
    
    // read the atoms
    for(int i=0; i < gtp.atoms.nr; i++) {
        gmx::t_atom *a;
        a = &(gtp.atoms.atom[i]);
        BeadType *type = top.GetOrCreateBeadType(*(gtp.atoms.atomtype[i]));
#ifdef GMX4CVS
        top.CreateBead(1, *(gtp.atoms.atomname[i]), type, a->resind, a->m, a->q);  
#else
        top.CreateBead(1, *(gtp.atoms.atomname[i]), type, a->resnr, a->m, a->q);  
#endif
        //cout << *(gtp.atoms.atomname[i]) << " residue: " << a->resnr << endl;
    }
    
    Bead *bead;
    for(int i=0; i<gtp.mols.nr; i++) {
        int i1 = gtp.mols.index[i];
        int i2 = gtp.mols.index[i+1];
        // if(i2-i1 < 2) continue;
        Molecule *mi = top.CreateMolecule("unnamed");    
        int res0 = top.getBead(i1)->getResnr();
        for(int i=i1; i<i2; ++i) {
            bead = top.getBead(i);
            stringstream n("");
            
            n << bead->getResnr()-res0 + 1 << ":" <<  top.getResidue(bead->getResnr())->getName() << ":" << bead->getName();
            //cout << n.str() << endl;
            mi->AddBead(top.getBead(i), n.str());
        }         
        
    }

    return true;
}
