/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
/// For earlier commit history see ctp commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#pragma once
#ifndef VOTCA_XTP_ATOM_H
#define VOTCA_XTP_ATOM_H

// Standard includes
#include "eigen.h"
#include <exception>
#include <map>
#include <string>
#include <votca/tools/types.h>
namespace votca {
namespace xtp {
class CptTable;
class Atom {
 public:
  // Real, direct MD-level bonded-partner atom IDs, up to this many --
  // covers virtually all ordinary organic chemistry (a typical carbon
  // has at most 4 bonded neighbors). A genuine atom with more bonded
  // partners than this (extremely rare) is a real, known limitation --
  // only the first kMaxBondedPartners found are recorded; this is not
  // guarded against/flagged with a hard error anywhere, deliberately,
  // to avoid making an already-rare edge case fatal.
  static constexpr Index kMaxBondedPartners = 4;

  struct data {
    Index id;
    char* element;
    char* name;
    double x;
    double y;
    double z;
    Index resnr;
    bool has_external_bond;
    double ext_bond_dir_x;
    double ext_bond_dir_y;
    double ext_bond_dir_z;
    Index ext_bond_partner_segment_id;
    Index bonded_partner_ids[kMaxBondedPartners];
  };
  Atom(Index resnr, std::string md_atom_name, Index atom_id,
       Eigen::Vector3d pos, std::string element);

  Atom(Index atom_id, std::string element, Eigen::Vector3d pos);

  Atom(data& d) { ReadData(d); }

  static std::string GetElementFromString(const std::string& MDName);

  Index getId() const { return id_; }
  const std::string& getName() const { return name_; }
  std::string getElement() const { return element_; }

  Index getResnr() const { return resnr_; }

  void setResnr(Index resnr) { resnr_ = resnr; }
  void Translate(const Eigen::Vector3d& shift) { pos_ = pos_ + shift; }

  void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& refPos);

  const Eigen::Vector3d& getPos() const { return pos_; }
  void setPos(const Eigen::Vector3d& r) { pos_ = r; }

  // Direction (normalized, unitless) toward an atom this one was
  // covalently bonded to at the MD level, but which fell outside the
  // mapped segment/fragment -- i.e. a bond that was "cut" by the
  // segment definition. Set directly by Md2QmEngine/SegmentMapper at
  // mapping time, using the exact same rigid-body transform already
  // applied to every other atom in the fragment (so this direction is
  // consistent with the fragment's own, already-idealized geometry,
  // not the raw, possibly out-of-plane/distorted MD-level direction).
  // A given atom may have at most one such recorded direction; an atom
  // with more than one external bond (rare) only retains the first one
  // set. Never set for atoms with no external bond at all (the normal
  // case for most atoms in most fragments) -- callers must check
  // hasExternalBond() before using getExternalBondDirection() at all.
  //
  // partner_atom_id is the raw, MD-level atom ID this external bond
  // points to (Atom::getId() of the partner) -- known and set here,
  // at the exact same point the direction itself is computed
  // (Md2QmEngine's own, real, MD-level bond connectivity). This is
  // deliberately transient/NOT persisted to the checkpoint file at
  // all (no CptTable column, unlike everything else on this class) --
  // it exists only to let Md2QmEngine look up, in a later, second
  // pass within the same map() call (once every segment's own ID is
  // known, which it is not yet at the point this direction is first
  // computed -- segments are still being built), which SEGMENT that
  // partner atom actually belongs to, and record that instead, via
  // setExternalBondPartnerSegmentId() below -- the actually useful,
  // persisted value. getExternalBondPartnerAtomId() is genuinely only
  // meant to be read within that same Md2QmEngine::map() call, never
  // afterward.
  bool hasExternalBond() const { return has_external_bond_; }
  const Eigen::Vector3d& getExternalBondDirection() const {
    return external_bond_direction_;
  }
  Index getExternalBondPartnerAtomId() const {
    return external_bond_partner_atom_id_;
  }
  void setExternalBondDirection(const Eigen::Vector3d& dir,
                                Index partner_atom_id) {
    has_external_bond_ = true;
    external_bond_direction_ = dir.normalized();
    external_bond_partner_atom_id_ = partner_atom_id;
  }

  // The SEGMENT ID (Segment::getId()) that this atom's own external-
  // bond partner (see above) belongs to -- -1 (this class's own
  // "-1 means unset" convention) until explicitly set. Unlike the
  // direction itself, a segment ID needs no rigid-body transform at
  // all (it is not a geometric quantity) -- SegmentMapper's own
  // TransferExternalBondDirection copies this straight across onto
  // the mapped QMAtom, alongside the direction transform, rather than
  // transforming it in any way. This IS persisted (a real CptTable
  // column, unlike getExternalBondPartnerAtomId() above) -- this is
  // the value later consumers (a planned linking-segment graph, and
  // the decision of whether a given external bond is already
  // satisfied within an assembled supermolecule) are actually meant
  // to read, from the checkpoint file, independently of
  // Md2QmEngine::map() itself having ever run in the same process.
  Index getExternalBondPartnerSegmentId() const {
    return external_bond_partner_segment_id_;
  }
  void setExternalBondPartnerSegmentId(Index segment_id) {
    external_bond_partner_segment_id_ = segment_id;
  }

  // Real, direct MD-level bonded-partner atom IDs (topological
  // connectivity, not geometry -- these are Atom::getId() values, MD
  // atom IDs, not array indices), for EVERY atom, not just ones with
  // an external bond -- unlike getExternalBondDirection() above, this
  // is meant to cover a fragment's own, full, internal connectivity
  // too (needed downstream for FragmentSaturator's own, planned
  // OpenBabel-based relaxation step, which needs real, known bond
  // connectivity for the whole fragment, not just the new bond, to
  // set up its own force field correctly at all -- see
  // FragmentSaturator's own header comment for why). Set directly by
  // Md2QmEngine at mapping time, from the real, MD-level topology's
  // own bond connectivity (the same source
  // getExternalBondDirection()'s own raw input comes from) -- these
  // IDs are NOT yet transformed/translated into QM-level IDs here;
  // that translation happens later, in SegmentMapper, matching the
  // same "raw here, transformed/translated later" split already used
  // for the external-bond direction. Unused slots hold -1, matching
  // this class's own, existing "-1 means unset" convention (id_'s own
  // default).
  const Index* getBondedPartnerIds() const { return bonded_partner_ids_; }
  void AddBondedPartner(Index partner_id) {
    for (Index& slot : bonded_partner_ids_) {
      if (slot == -1) {
        slot = partner_id;
        return;
      }
    }
    // All kMaxBondedPartners slots already full -- see this class's
    // own, static kMaxBondedPartners comment above for why this is
    // silently dropped rather than treated as a hard error.
  }

  std::string identify() const { return "atom"; }

  friend std::ostream& operator<<(std::ostream& out, const Atom& atom) {
    out << atom.getId() << " " << atom.getName() << " " << atom.getElement()
        << " " << atom.getResnr();
    out << " " << atom.getPos().x() << "," << atom.getPos().y() << ","
        << atom.getPos().z() << "\n";
    return out;
  }

  static void SetupCptTable(CptTable& table);

  void WriteData(data& d) const;

  void ReadData(const data& d);

 private:
  Index id_ = -1;
  std::string name_ = "";

  std::string element_ = "";
  Index resnr_ = -1;
  Eigen::Vector3d pos_ = Eigen::Vector3d::Zero();
  bool has_external_bond_ = false;
  Eigen::Vector3d external_bond_direction_ = Eigen::Vector3d::Zero();
  Index external_bond_partner_atom_id_ = -1;
  Index external_bond_partner_segment_id_ = -1;
  Index bonded_partner_ids_[kMaxBondedPartners] = {-1, -1, -1, -1};
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ATOM_H
