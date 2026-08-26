/*
 *            Copyright 2009-2023 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_QMMOLECULE_H
#define VOTCA_XTP_QMMOLECULE_H

// Local VOTCA includes
#include "atomcontainer.h"
#include "votca/xtp/ewald/polarseg.h"
#include "qmatom.h"

namespace votca {
namespace xtp {

class QMMolecule : public AtomContainer<QMAtom> {
 public:
  QMMolecule(std::string name, Index id) : AtomContainer<QMAtom>(name, id) {};

  QMMolecule(CheckpointReader& r) : AtomContainer<QMAtom>(r) {};
  void LoadFromFile(std::string filename);

  void FromPolarSegment( PolarSeg* polar);

  void WriteXYZ(std::string filename, std::string header) const;

  void AddContainer(const AtomContainer<QMAtom>& container) {
    Index offset = atomlist_.size();
    type_ += "_" + container.getType();
    for (const auto& at : container) {
      // Real, direct fix for a real, direct, pre-existing bug --
      // reconstructing a brand new QMAtom from scratch here (the
      // (Index, element, pos) constructor, as this used to do)
      // silently discards everything else about the original atom
      // at all: hasExternalBond()/external_bond_direction_/
      // external_bond_partner_segment_id_/bonded_partner_ids_ all
      // reset to their own, empty defaults. Confirmed directly, this
      // exact way, from the user's own real, direct diagnostic run:
      // only ONE of three genuinely expected external bonds ever
      // actually saturated at all -- specifically the one on the
      // atom belonging to the FIRST segment merged into a QMMolecule
      // this way (qmmol = mapper.map(*seg1, stateA), a real, direct
      // copy-construction, not affected by this bug at all) -- every
      // subsequent segment's own atoms, merged in via THIS method
      // (qmmol.AddContainer(mapper.map(seg2, stateB))), silently lost
      // this data entirely.
      //
      // Fixed by copying the whole atom directly (preserving
      // everything about it), then only mutating what genuinely does
      // need to change once merged into a larger container: its own
      // id (via setID(), matching this method's own, already-
      // established "unique ids" comment/purpose exactly), and its
      // own bonded_partner_ids_ (via the new setBondedPartnerIds()
      // added directly alongside this fix, qmatom.h) -- these are
      // recorded local to the smaller QMMolecule this atom was
      // originally mapped within, so need the exact same offset
      // applied to them too, or they would end up pointing at the
      // wrong atom entirely once merged. hasExternalBond()/
      // external_bond_direction_/external_bond_partner_segment_id_
      // need no such offsetting at all -- none of them are ids that
      // point at another atom within this same container (the first
      // two are per-atom properties; the third is a real, direct
      // SEGMENT id, an entirely separate, global id space).
      QMAtom atom = at;
      atom.setID(at.getId() + offset);
      const Index* partners = at.getBondedPartnerIds();
      Index offset_partners[QMAtom::kMaxBondedPartners];
      for (Index i = 0; i < QMAtom::kMaxBondedPartners; i++) {
        offset_partners[i] = (partners[i] == -1) ? -1 : partners[i] + offset;
      }
      atom.setBondedPartnerIds(offset_partners);
      atomlist_.push_back(atom);
    }
    calcPos();
  }

  void ReorderAtomIDs() {
    Index id = 0;
    for (auto& at : atomlist_) {
      at.setID(id);
      id++;
    }
  }

  friend std::ostream& operator<<(std::ostream& out,
                                  const QMMolecule& container) {
    out << container.getId() << " " << container.getType() << "\n";
    for (const QMAtom& atom : container) {
      out << atom;
    }
    out << std::endl;
    return out;
  }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMMOLECULE_H
