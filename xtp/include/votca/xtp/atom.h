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
  bool hasExternalBond() const { return has_external_bond_; }
  const Eigen::Vector3d& getExternalBondDirection() const {
    return external_bond_direction_;
  }
  void setExternalBondDirection(const Eigen::Vector3d& dir) {
    has_external_bond_ = true;
    external_bond_direction_ = dir.normalized();
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
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ATOM_H
