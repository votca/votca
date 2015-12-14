/* 
 * Copyright 2009-2015 The VOTCA Development Team (http://www.votca.org)
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

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include <iostream>
#include "gmxtopologyreader.h"

#if (GMX == 51)||(GMX == 52)
        #include <gromacs/fileio/tpxio.h>
        #include <gromacs/topology/atoms.h>
        #include <gromacs/topology/topology.h>
#elif GMX == 50
        #include <gromacs/fileio/tpxio.h>
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

    t_inputrec ir;
    ::matrix gbox;

#if GMX == 52
    (void)read_tpx((char *)file.c_str(),&ir,gbox,&natoms,NULL,NULL,&mtop);
#else
    (void)read_tpx((char *)file.c_str(),&ir,gbox,&natoms,NULL,NULL,NULL,&mtop);
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

    int ifirstatom = 0;
    for(int iblock=0; iblock<mtop.nmolblock; ++iblock) {
        gmx_moltype_t *mol
                = &(mtop.moltype[mtop.molblock[iblock].type]);

        string molname =  *(mol->name);

        int res_offset = top.ResidueCount();

        t_atoms *atoms=&(mol->atoms);

        for(int i=0; i < atoms->nres; i++) {
                top.CreateResidue(*(atoms->resinfo[i].name));
        }

        for(int imol=0; imol<mtop.molblock[iblock].nmol; ++imol) {
            Molecule *mi = top.CreateMolecule(molname);

            // read the atoms
            for(int iatom=0; iatom<mtop.molblock[iblock].natoms_mol; iatom++) {
                t_atom *a = &(atoms->atom[iatom]);

                BeadType *type = top.GetOrCreateBeadType(*(atoms->atomtype[iatom]));
                Bead *bead = top.CreateBead(1, *(atoms->atomname[iatom]), type, a->resind + res_offset, a->m, a->q);

                stringstream nm;
                nm << bead->getResnr() + 1 - res_offset << ":" <<  top.getResidue(bead->getResnr())->getName() << ":" << bead->getName();
                mi->AddBead(bead, nm.str());
            }

            // add exclusions
            for(int iatom=0; iatom<mtop.molblock[iblock].natoms_mol; iatom++) {
                // read exclusions
                t_blocka * excl = &(mol->excls);
                // insert exclusions
                list<Bead *> excl_list;
                for(int k=excl->index[iatom]; k<excl->index[iatom+1]; k++) {
                	excl_list.push_back(top.getBead(excl->a[k]+ifirstatom));
                }
                top.InsertExclusion(top.getBead(iatom+ifirstatom), excl_list);
            }
            ifirstatom+=mtop.molblock[iblock].natoms_mol;
        }
    }

    matrix m;
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            m[i][j] = gbox[j][i];
    top.setBox(m);

    return true;
}

}}
