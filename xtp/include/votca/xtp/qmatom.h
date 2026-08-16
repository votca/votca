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
  Index bonded_partner_ids_[kMaxBondedPartners] = {-1, -1, -1, -1};

 public:
  static void SetupCptTable(CptTable& table);

  void WriteData(data& d) const;

  void ReadData(const data& d);
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMATOM_H
