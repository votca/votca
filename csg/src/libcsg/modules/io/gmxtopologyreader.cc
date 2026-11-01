/*
 * Copyright 2009-2023 The VOTCA Development Team (http://www.votca.org)
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

#include <array>
#include <iostream>
#include <string>

#include "../../../../include/votca/csg/interaction.h"
#include "../../../../include/votca/csg/topology.h"

#include "gmxtopologyreader.h"

#include <gromacs/fileio/tpxio.h>
#include <gromacs/mdtypes/inputrec.h>
#include <gromacs/topology/atoms.h>
#include <gromacs/topology/ifunc.h>
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

  size_t nmolblock = mtop.molblock.size();

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

      size_t natoms_mol = mtop.moltype[mtop.molblock[iblock].type].atoms.nr;
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
           << top.getResidue(bead->getResnr()).getName() << ":"
           << bead->getName();
        mi->AddBead(bead, nm.str());
      }

      // add exclusions
      for (size_t iatom = 0; iatom < natoms_mol; iatom++) {
        std::list<Bead *> excl_list;
        gmx::ListOfLists<int> &excl = mol->excls;
        for (const Index k : excl[iatom]) {
          excl_list.push_back(top.getBead(k + ifirstatom));
        }
        top.InsertExclusion(top.getBead(iatom + ifirstatom), excl_list);
      }

      // GROMACS 2026 and later only provide the enum class
      // InteractionFunction; 2023-2025 only provide the plain enum values
      // F_BONDS/F_CONSTR. GMX_VERSION is from gromacs/version.h.
#if GMX_VERSION >= 20260000
      std::array<int, 2> ftypes = {
          static_cast<int>(InteractionFunction::Bonds),
          static_cast<int>(InteractionFunction::Constraints)};
#else
      std::array<int, 2> ftypes = {F_BONDS, F_CONSTR};
#endif

      // mol->ilist holds one flat list per interaction type; each entry
      // occupies interaction_function[ftype].nratoms + 1 ints: the
      // parameter type index followed by the atom indices. Bonds and
      // constraints are read, since "constraints = h-bonds/all-bonds"
      // turns bonds into constraints. Other bond-like types (G96, Morse,
      // cubic bonds, connections) are not read.
      for (int ftype : ftypes) {
        const InteractionList &ilist = mol->ilist[ftype];
        Index nratoms = interaction_function[ftype].nratoms;
        Index stride = nratoms + 1;
        for (Index i = 0; i < Index(ilist.size()); i += stride) {
          // ilist.iatoms[i] is the parameter index (not needed, only
          // connectivity is read).
          Index atom1 = ilist.iatoms[i + 1];
          Index atom2 = ilist.iatoms[i + 2];
          // getGroup() asserts a non-empty group, which a new IBond does
          // not have; use "BONDS" as in lammpsdatareader.cc and
          // pdbreader.cc.
          Interaction *ic =
              new IBond(Index(atom1 + ifirstatom), Index(atom2 + ifirstatom));
          ic->setGroup("BONDS");
          top.AddBondedInteraction(ic);
        }
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
