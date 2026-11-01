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
  // Maximum number of stored bonded partners (sufficient for organic
  // molecules). Further partners are silently dropped rather than treated
  // as an error.
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

  // Unit vector toward a covalently bonded MD atom that lies outside the
  // segment, i.e. a bond cut by the segment definition. Md2QmEngine sets it
  // from the MD geometry; SegmentMapper rotates it with the fragment's
  // MD->template transform for the mapped QMAtom. Only the first external
  // bond of an atom is stored; check hasExternalBond() before use.
  //
  // partner_atom_id is the MD atom ID of the partner. It is transient (not
  // written to the checkpoint file) and only used within
  // Md2QmEngine::map() to determine the partner's segment ID
  // (setExternalBondPartnerSegmentId()) once all segments exist.
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

  // Segment::getId() of the segment the external-bond partner belongs to,
  // -1 if unset. Copied unchanged to the mapped QMAtom by SegmentMapper.
  // Persisted in the checkpoint file, e.g. for the segment bond graph used
  // by Topology::FindLinkingSegments.
  Index getExternalBondPartnerSegmentId() const {
    return external_bond_partner_segment_id_;
  }
  void setExternalBondPartnerSegmentId(Index segment_id) {
    external_bond_partner_segment_id_ = segment_id;
  }

  // IDs (Atom::getId(), not indices) of all covalently bonded atoms, set by
  // Md2QmEngine from the MD topology for every atom. The fragment relaxation
  // in FragmentSaturator needs this full internal connectivity. The IDs are
  // MD-level here and translated to QM-level IDs in SegmentMapper. Unused
  // slots hold -1.
  const Index* getBondedPartnerIds() const { return bonded_partner_ids_; }
  void AddBondedPartner(Index partner_id) {
    for (Index& slot : bonded_partner_ids_) {
      if (slot == -1) {
        slot = partner_id;
        return;
      }
    }
    // All slots full: partner is dropped (see kMaxBondedPartners).
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
