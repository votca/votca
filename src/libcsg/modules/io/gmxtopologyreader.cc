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
    gmx::gmx_mtop_t mtop;
    char       title[STRLEN];
    gmx::rvec       *xtop;
    gmx::matrix     box;

    // whateverer are these
    int sss;
    int natoms;
    gmx::real    ttt,lll;
    // cleanup topology to store new data
    top.Cleanup();

    int ePBC =
        read_tpx((char *)file.c_str(),&sss,&ttt,&lll,NULL,NULL,&natoms,NULL,NULL,NULL,&mtop);

    int count=0;
    for(int iblock=0; iblock<mtop.nmolblock; ++iblock)
        count+=mtop.molblock[iblock].nmol;

    if(count != mtop.mols.nr  ) {
        throw runtime_error("gromacs topology contains inconsistency in molecule definitons\n\n"
                "A possible reason is an outdated .tpr file. Please rerun grompp to generate a new tpr file.\n"
                "If the problem remains or "
                "you're missing the files to rerun grompp,\ncontact the votca mailing list for a solution.");
    }

    for(int iblock=0; iblock<mtop.nmolblock; ++iblock) {
        gmx::gmx_moltype_t *mol
                = &(mtop.moltype[mtop.molblock[iblock].type]);

        string molname =  *(mol->name);

        int res_offset = top.ResidueCount();

        gmx::t_atoms *atoms=&(mol->atoms);

        for(int i=0; i < atoms->nres; i++) {
//            #ifdef GMX4CVS
//                top.CreateResidue(*(atoms->resinfo[i].name));
//            #else
                top.CreateResidue(*(atoms->resname[i]));
//            #endif
        }


        for(int imol=0; imol<mtop.molblock[iblock].nmol; ++imol) {
            Molecule *mi = top.CreateMolecule(molname);

            // read the atoms
            for(int iatom=0; iatom<mtop.molblock[iblock].natoms_mol; iatom++) {
                gmx::t_atom *a = &(atoms->atom[iatom]);

                BeadType *type = top.GetOrCreateBeadType(*(atoms->atomtype[iatom]));
//            #ifdef GMX4CVS
//                Bead *bead = top.CreateBead(1, *(atoms->atomname[iatom]), type, a->resind, a->m, a->q);
//            #else
                Bead *bead = top.CreateBead(1, *(atoms->atomname[iatom]), type, a->resnr, a->m, a->q);
//            #endif

                stringstream nm;
                nm << bead->getResnr() + 1 << ":" <<  top.getResidue(res_offset + bead->getResnr())->getName() << ":" << bead->getName();
                mi->AddBead(bead, nm.str());
            }
        }
    }

    return true;
}
