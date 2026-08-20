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
#ifndef VOTCA_XTP_QMATOM_H
#define VOTCA_XTP_QMATOM_H

// VOTCA includes
#include "eigen.h"
#include <votca/tools/elements.h>
#include <votca/tools/types.h>

namespace votca {
namespace xtp {
class CptTable;
/**
 *    \brief container for QM atoms
 *
 *    Stores atom type, coordinates, charge
 */
class QMAtom {
  friend class ECPAOBasis;

 public:
  // Same value as Atom::kMaxBondedPartners -- see its own comment
  // (atom.h) for why 4.
  static constexpr Index kMaxBondedPartners = 4;

  struct data {
    Index index;
    char* element;
    double x;
    double y;
    double z;
    Index nuccharge;
    Index ecpcharge;
    bool has_external_bond;
    double ext_bond_dir_x;
    double ext_bond_dir_y;
    double ext_bond_dir_z;
    Index ext_bond_partner_segment_id;
    Index bonded_partner_ids[kMaxBondedPartners];
  };

  QMAtom(Index index, std::string element, Eigen::Vector3d pos);

  QMAtom(const data& d);

  const Eigen::Vector3d& getPos() const { return pos_; }

  void Translate(const Eigen::Vector3d& shift) { pos_ += shift; }

  void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& refPos);

  void setPos(const Eigen::Vector3d& position) { pos_ = position; }

  const std::string& getElement() const { return element_; }

  Index getId() const { return index_; }

  void setID(const Index index) { index_ = index; }

  Index getNuccharge() const { return nuccharge_ - ecpcharge_; }

  Index getElementNumber() const { return nuccharge_; }

  // Same meaning as Atom::hasExternalBond()/getExternalBondDirection()
  // (a direction toward an MD-level atom this one was covalently
  // bonded to, but which fell outside the mapped fragment) -- a
  // SEPARATE field from Atom's own, since QMAtom does not inherit
  // from Atom at all. Set directly by SegmentMapper::MapMapAtomonMD/
  // PlaceMapAtomonMD, by transforming the corresponding MD-level
  // Atom's own, already-recorded, raw direction through the exact
  // same rigid-body transform already applied to every other atom in
  // the fragment there -- never computed independently here.
  bool hasExternalBond() const { return has_external_bond_; }
  const Eigen::Vector3d& getExternalBondDirection() const {
    return external_bond_direction_;
  }
  void setExternalBondDirection(const Eigen::Vector3d& dir) {
    has_external_bond_ = true;
    external_bond_direction_ = dir.normalized();
  }

  // Marks this atom's own external bond as resolved/no longer
  // needing saturation -- e.g. because, within a specific, assembled
  // supermolecule, its own external-bond partner segment
  // (getExternalBondPartnerSegmentId()) turns out to already be
  // present too, so the bond is already satisfied there and does not
  // need a new H at all. Deliberately resets has_external_bond_ to
  // false entirely (not merely a separate "skip saturation" flag) --
  // once a bond is known to be satisfied within a given, specific
  // supermolecule, it genuinely is no longer "external" to it at all,
  // so hasExternalBond() itself should honestly report false for it,
  // not some separate, parallel notion of "external but exempt".
  // Deliberately only ever affects THIS, specific, already-mapped
  // QMAtom instance/copy -- SegmentMapper::map() constructs a
  // brand-new QMMolecule from scratch on every call, re-deriving
  // hasExternalBond()/getExternalBondPartnerSegmentId() fresh each
  // time from the underlying MD-level Segment's own, persisted atom
  // data (never from a previously-cleared copy) -- so clearing this
  // on one, specific mapped copy (e.g. within one specific pair's own
  // supermolecule) has no effect at all on any other, separate
  // mapping of the same underlying segment (e.g. within a different
  // pair).
  void clearExternalBond() {
    has_external_bond_ = false;
    external_bond_direction_ = Eigen::Vector3d::Zero();
    external_bond_partner_segment_id_ = -1;
  }

  // Same meaning as Atom::getExternalBondPartnerSegmentId() -- the
  // Segment::getId() this external bond crosses into. Unlike the
  // direction itself, a segment id needs no rigid-body transform at
  // all (it is not a geometric quantity) -- SegmentMapper copies it
  // straight across from the corresponding MD-level Atom, alongside
  // the direction transform, rather than transforming it in any way.
  // Only meaningful when hasExternalBond() is true (matches -1, this
  // class's own "-1 means unset" convention, otherwise). Unlike
  // Atom's own version, QMAtom needs no separate, transient "partner
  // ATOM id" field at all -- that value only ever existed to let
  // Md2QmEngine resolve it into this, real, final segment id in the
  // first place (see Atom::getExternalBondPartnerAtomId's own header
  // comment, atom.h); QMAtom only ever needs the already-resolved
  // result.
  Index getExternalBondPartnerSegmentId() const {
    return external_bond_partner_segment_id_;
  }
  void setExternalBondPartnerSegmentId(Index segment_id) {
    external_bond_partner_segment_id_ = segment_id;
  }

  // Same meaning as Atom::getBondedPartnerIds() (a fragment's own,
  // full, internal bond connectivity) -- but, unlike that one, these
  // are already QM-LEVEL atom IDs (i.e. QMAtom::getId() values,
  // matching this same QMMolecule's own atom numbering), NOT the raw,
  // MD-level IDs Atom itself stores. Set directly by
  // SegmentMapper::MapMapAtomonMD/PlaceMapAtomonMD, by translating
  // each corresponding MD-level Atom's own, already-recorded, raw
  // partner IDs into QM-level ones there, via the same
  // mapatom_ids/mdatom_ids correspondence already used elsewhere in
  // that same class -- never computed independently here. Same
  // "-1 means unset" sentinel convention as Atom's own.
  const Index* getBondedPartnerIds() const { return bonded_partner_ids_; }
  void AddBondedPartner(Index partner_id) {
    for (Index& slot : bonded_partner_ids_) {
      if (slot == -1) {
        slot = partner_id;
        return;
      }
    }
    // All kMaxBondedPartners slots already full -- see
    // Atom::kMaxBondedPartners's own comment for why this is silently
    // dropped rather than treated as a hard error.
  }
  // Real, direct addition -- needed for QMMolecule::AddContainer()
  // (qmmolecule.h) to be able to correctly offset an atom's own
  // bonded_partner_ids_ when merging it into a larger QMMolecule
  // (its own ids, as recorded, are local to the smaller QMMolecule
  // it was originally mapped within -- a real, direct offset is
  // needed once merged into a bigger one). AddBondedPartner() alone
  // cannot do this: it only ever appends into the next available -1
  // slot, never overwrites the whole array at once, so it cannot
  // replace already-set values with their own, new, offset ones.
  void setBondedPartnerIds(const Index (&ids)[kMaxBondedPartners]) {
    for (Index i = 0; i < kMaxBondedPartners; i++) {
      bonded_partner_ids_[i] = ids[i];
    }
  }

  std::string identify() const { return "qmatom"; }

  friend std::ostream& operator<<(std::ostream& out, const QMAtom& atom) {
    out << atom.getId() << " " << atom.getElement();
    out << " " << atom.getPos().x() << "," << atom.getPos().y() << ","
        << atom.getPos().z() << " " << atom.getNuccharge() << "\n";
    return out;
  }

 private:
  Index index_;
  std::string element_;
  Eigen::Vector3d pos_;  // Bohr
  Index nuccharge_ = 0;
  Index ecpcharge_ = 0;  // ecp charge is set in ecpaobasis.fill
  bool has_external_bond_ = false;
  Eigen::Vector3d external_bond_direction_ = Eigen::Vector3d::Zero();
  Index external_bond_partner_segment_id_ = -1;
  Index bonded_partner_ids_[kMaxBondedPartners] = {-1, -1, -1, -1};

 public:
  static void SetupCptTable(CptTable& table);

  void WriteData(data& d) const;

  void ReadData(const data& d);
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMATOM_H
