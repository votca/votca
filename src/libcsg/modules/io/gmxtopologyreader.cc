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
#include "gmxtopologyreader.h"
#include "votca_config.h"

#if GMX == 50
        #include <gromacs/legacyheaders/statutil.h>
        #include <gromacs/legacyheaders/typedefs.h>
        #include <gromacs/legacyheaders/smalloc.h>
        #include <gromacs/legacyheaders/vec.h>
        #include <gromacs/legacyheaders/copyrite.h>
        #include <gromacs/legacyheaders/statutil.h>
        #include <gromacs/legacyheaders/tpxio.h>
#elif GMX == 45
        #include <gromacs/statutil.h>
        #include <gromacs/typedefs.h>
        #include <gromacs/smalloc.h>
        #include <gromacs/vec.h>
        #include <gromacs/copyrite.h>
        #include <gromacs/statutil.h>
        #include <gromacs/tpxio.h>
#elif GMX == 40
    extern "C"
    {
        #include <statutil.h>
        #include <typedefs.h>
        #include <smalloc.h>
        #include <vec.h>
        #include <copyrite.h>
        #include <statutil.h>
        #include <tpxio.h>
    }
#else
#error Unsupported GMX version
#endif
    // this one is needed because of bool is defined in one of the headers included by gmx
    #undef bool

namespace votca { namespace csg {

bool GMXTopologyReader::ReadTopology(string file, Topology &top)
{ 
    gmx_mtop_t mtop;

    int natoms;
    // cleanup topology to store new data
    top.Cleanup();
    set_program_name("VOTCA");

#if GMX == 50
    t_inputrec ir;
    ::matrix gbox;

    (void)read_tpx((char *)file.c_str(),&ir,gbox,&natoms,NULL,NULL,NULL,&mtop);
#elif GMX == 45
    t_inputrec ir;
    ::matrix gbox;
    (void)read_tpx((char *)file.c_str(),&ir,gbox,&natoms,NULL,NULL,NULL,&mtop);
#elif GMX == 40
    int sss;   // wtf is this
    ::real    ttt,lll; // wtf is this
    (void)read_tpx((char *)file.c_str(),&sss,&ttt,&lll,NULL,NULL,&natoms,NULL,NULL,NULL,&mtop);
#else
#error Unsupported GMX version
#endif

    int count=0;
    for(int iblock=0; iblock<mtop.nmolblock; ++iblock)
        count+=mtop.molblock[iblock].nmol;

    if(count != mtop.mols.nr  ) {
        throw runtime_error("gromacs topology contains inconsistency in molecule definitons\n\n"
                "A possible reason is an outdated .tpr file. Please rerun grompp to generate a new tpr file.\n"
                "If the problem remains or "
                "you're missing the files to rerun grompp,\n contact the votca mailing list for a solution.");
    }

    for(int iblock=0; iblock<mtop.nmolblock; ++iblock) {
        gmx_moltype_t *mol
                = &(mtop.moltype[mtop.molblock[iblock].type]);

        string molname =  *(mol->name);

        int res_offset = top.ResidueCount();

        t_atoms *atoms=&(mol->atoms);

        for(int i=0; i < atoms->nres; i++) {
#if GMX == 50
                top.CreateResidue(*(atoms->resinfo[i].name));
#elif GMX == 45
                top.CreateResidue(*(atoms->resinfo[i].name));
#elif GMX == 40
                top.CreateResidue(*(atoms->resname[i]));
#else
#error Unsupported GMX version
#endif
        }

        int ifirstatom = 0;
        for(int imol=0; imol<mtop.molblock[iblock].nmol; ++imol) {
            Molecule *mi = top.CreateMolecule(molname);

            // read the atoms
            for(int iatom=0; iatom<mtop.molblock[iblock].natoms_mol; iatom++) {
                t_atom *a = &(atoms->atom[iatom]);

                // read exclusions
                t_blocka * excl = &(mol->excls);
                // insert exclusions
                list<int> excl_list;
                for(int k=excl->index[iatom]; k<excl->index[iatom+1]; k++) {
                    excl_list.push_back(excl->a[k]+ifirstatom);
                }
                top.InsertExclusion(iatom, excl_list);

                BeadType *type = top.GetOrCreateBeadType(*(atoms->atomtype[iatom]));
#if GMX == 50
                Bead *bead = top.CreateBead(1, *(atoms->atomname[iatom]), type, a->resind, a->m, a->q);
#elif GMX == 45
                Bead *bead = top.CreateBead(1, *(atoms->atomname[iatom]), type, a->resind, a->m, a->q);
#elif GMX == 40
                Bead *bead = top.CreateBead(1, *(atoms->atomname[iatom]), type, a->resnr, a->m, a->q);
#else
#error Unsupported GMX version
#endif

                stringstream nm;
                nm << bead->getResnr() + 1 << ":" <<  top.getResidue(res_offset + bead->getResnr())->getName() << ":" << bead->getName();
                mi->AddBead(bead, nm.str());
            }
            ifirstatom+=mtop.molblock[iblock].natoms_mol;
        }
    }

    return true;
}

}}
