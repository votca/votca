/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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
#include <string>

#include "../../../../include/votca/csg/topology.h"

#include "gmxtopologyreader.h"

#include <gromacs/fileio/tpxio.h>
#include <gromacs/mdtypes/inputrec.h>
#include <gromacs/topology/atoms.h>
#include <gromacs/topology/topology.h>
#include <gromacs/version.h>

// this one is needed because of bool is defined in one of the headers included
// by gmx
#undef bool

namespace votca {
namespace csg {

bool GMXTopologyReader::ReadTopology(std::string file, Topology &top) {
  gmx_mtop_t mtop;

  int natoms;
  // cleanup topology to store new data
  top.Cleanup();

  t_inputrec ir;
  ::matrix gbox;

  (void)read_tpx((char *)file.c_str(), &ir, gbox, &natoms, nullptr, nullptr,
                 &mtop);

  size_t ifirstatom = 0;

#if GMX_VERSION >= 20190000
  size_t nmolblock = mtop.molblock.size();
#else
  size_t nmolblock = mtop.nmolblock;
#endif

  for (size_t iblock = 0; iblock < nmolblock; ++iblock) {
    gmx_moltype_t *mol = &(mtop.moltype[mtop.molblock[iblock].type]);

    std::string molname = *(mol->name);

    Index res_offset = top.ResidueCount();

    t_atoms *atoms = &(mol->atoms);

    for (Index i = 0; i < atoms->nres; i++) {
      top.CreateResidue(*(atoms->resinfo[i].name));
    }

    for (Index imol = 0; imol < mtop.molblock[iblock].nmol; ++imol) {
      Molecule *mi = top.CreateMolecule(molname);

#if GMX_VERSION >= 20190000
      size_t natoms_mol = mtop.moltype[mtop.molblock[iblock].type].atoms.nr;
#else
      size_t natoms_mol = mtop.molblock[iblock].natoms_mol;
#endif
      // read the atoms
      for (size_t iatom = 0; iatom < natoms_mol; iatom++) {
        t_atom *a = &(atoms->atom[iatom]);

        std::string bead_type = *(atoms->atomtype[iatom]);
        if (!top.BeadTypeExist(bead_type)) {
          top.RegisterBeadType(bead_type);
        }
        Bead *bead =
            top.CreateBead(Bead::spherical, *(atoms->atomname[iatom]),
                           bead_type, a->resind + res_offset, a->m, a->q);

        std::stringstream nm;
        nm << bead->getResnr() + 1 - res_offset << ":"
           << top.getResidue(bead->getResnr())->getName() << ":"
           << bead->getName();
        mi->AddBead(bead, nm.str());
      }

      // add exclusions
      for (size_t iatom = 0; iatom < natoms_mol; iatom++) {
        // read exclusions
        t_blocka *excl = &(mol->excls);
        // insert exclusions
        std::list<Bead *> excl_list;
        for (Index k = excl->index[iatom]; k < excl->index[iatom + 1]; k++) {
          excl_list.push_back(top.getBead(excl->a[k] + ifirstatom));
        }
        top.InsertExclusion(top.getBead(iatom + ifirstatom), excl_list);
      }
      ifirstatom += natoms_mol;
    }
  }

  Eigen::Matrix3d m;
  for (Index i = 0; i < 3; i++) {
    for (Index j = 0; j < 3; j++) {
      m(i, j) = gbox[j][i];
    }
  }
  top.setBox(m);

  return true;
}

}  // namespace csg
}  // namespace votca
