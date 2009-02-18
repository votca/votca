// 
// File:   pdbtopologyreader.cc
// Author: victor
//
// Created on 27. Dezember 2007, 21:19
//

#include <iostream>
#include "pdbtopologyreader.h"

namespace gmx {
   extern "C" {
        #include <statutil.h>
        #include <typedefs.h>
        #include <smalloc.h>
        #include <confio.h>
        #include <vec.h>
        #include <copyrite.h>
        #include <statutil.h>
        #include <tpxio.h>
    }
    // this one is needed because of bool is defined in one of the headers included by gmx
    #undef bool
}

#define snew2(ptr,nelem) (ptr)=(gmx::rvec*)gmx::save_calloc(#ptr,__FILE__,__LINE__,\
                        (nelem),sizeof(*(ptr)))


bool PDBTopologyReader::ReadTopology(string file, Topology &top)
{
    char title[512];
    gmx::rvec *x, *v;
    gmx::matrix box;
    int ePBC;
    gmx::t_atoms atoms;
    
    //snew(atoms,1);
    gmx::get_stx_coordnum((char*)file.c_str(),&(atoms.nr));
    init_t_atoms(&atoms,atoms.nr,TRUE);
    snew2(x,atoms.nr);
    snew2(v,atoms.nr);
    fprintf(stderr,"\nReading structure file\n");

    gmx::read_stx_conf((char*)file.c_str(), title,&atoms,
                          x,NULL,&ePBC,box);
    Residue *res = top.CreateResidue("no");
    // read the atoms
    for(int i=0; i < atoms.nr; i++) {
        gmx::t_atom *a;
        a = &(atoms.atom[i]);
        BeadType *type = top.GetOrCreateBeadType("no");
        top.CreateBead(1, *(atoms.atomname[i]), type, res->getId(), a->m, a->q);  
        //cout << *(gtp.atoms.atomname[i]) << " residue: " << a->resnr << endl;
    }
   
    return true;
}
