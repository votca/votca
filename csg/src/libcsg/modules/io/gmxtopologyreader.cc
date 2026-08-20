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

      // Real, direct, root-cause fix -- this function used to never
      // read any real bond connectivity at all, from any .tpr file,
      // regardless of its own actual content (confirmed directly,
      // this same session, by reading this function in full, before
      // this fix, and finding no mention of bonds/interactions
      // anywhere in it at all) -- meaning top.BondedInteractions()
      // was always empty for any topology loaded this way, silently
      // breaking every downstream consumer relying on real bond
      // connectivity (e.g. xtp's own external-bond-direction
      // detection, Md2QmEngine::map, md2qmengine.cc, built and
      // extensively tested earlier, against a hand-constructed bond
      // map -- never against a real, actual GMXTopologyReader-loaded
      // topology at all, until the user's own real run surfaced this
      // real, direct gap).
      //
      // mol->ilist (InteractionLists, a real GROMACS type) holds one
      // real, flat InteractionList per interaction type, indexed by
      // InteractionFunction -- each real, individual interaction
      // occupies interaction_function[ft].nratoms + 1 consecutive
      // int entries within its own list's own iatoms array: the
      // interaction TYPE index first, then nratoms real atom
      // indices. Only InteractionFunction::Bonds and
      // InteractionFunction::Constraints are read here -- covers
      // ordinary bonds (regardless of any .mdp constraints setting;
      // "constraints = h-bonds"/"all-bonds" converts some real bonds
      // into Constraints entries specifically instead of Bonds, so
      // both are checked) -- but NOT every other, less common,
      // force-field-specific "bond-like" interaction type GROMACS
      // itself supports (e.g. G96Bonds, Morse, Cubic, Connections) --
      // a real, direct, honest limitation, worth being aware of for
      // force fields that use one of those instead of ordinary
      // harmonic bonds.
      for (InteractionFunction ft :
          {InteractionFunction::Bonds, InteractionFunction::Constraints}) {
        const InteractionList &ilist = mol->ilist[ft];
        Index nratoms = interaction_function[static_cast<int>(ft)].nratoms;
        Index stride = nratoms + 1;
        for (Index i = 0; i < Index(ilist.size()); i += stride) {
          // ilist.iatoms[i] itself is the interaction TYPE index
          // (into mtop.ffparams.iparams) -- not an atom index at
          // all, and not needed here, since only real connectivity
          // (which atoms are bonded), not the real bond's own force
          // constant/equilibrium length, is wanted at this level.
          Index atom1 = ilist.iatoms[i + 1];
          Index atom2 = ilist.iatoms[i + 2];
          top.AddBondedInteraction(
              new IBond(Index(atom1 + ifirstatom), Index(atom2 + ifirstatom)));
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
