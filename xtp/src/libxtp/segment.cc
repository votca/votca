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

// Local VOTCA includes
#include "votca/xtp/segment.h"
#include "votca/xtp/atomcontainer.h"

namespace votca {
namespace xtp {

double Segment::getApproxSize() const {
  std::pair<Eigen::Vector3d, Eigen::Vector3d> minmax = CalcSpatialMinMax();
  return (minmax.first - minmax.second).norm();
}

const Atom* Segment::getAtom(Index id) const {

  for (const Atom& atom : *this) {
    if (atom.getId() == id) {
      return &atom;
    }
  }
  return nullptr;
}

void Segment::WriteToCpt(CheckpointWriter& w) const {
  AtomContainer<Atom>::WriteToCpt(w);
  w(molecule_ids_, "mol_ids");

  w(U_xX_nN_.getValue(QMStateType::Electron), "U_xX_nN_e");
  w(U_xX_nN_.getValue(QMStateType::Hole), "U_xX_nN_h");
  w(U_xX_nN_.getValue(QMStateType::Singlet), "U_xX_nN_s");
  w(U_xX_nN_.getValue(QMStateType::Triplet), "U_xX_nN_t");
  w(U_nX_nN_.getValue(QMStateType::Electron), "U_nX_nN_e");
  w(U_nX_nN_.getValue(QMStateType::Hole), "U_nX_nN_h");
  w(U_nX_nN_.getValue(QMStateType::Singlet), "U_nX_nN_s");
  w(U_nX_nN_.getValue(QMStateType::Triplet), "U_nX_nN_t");
  w(U_xN_xX_.getValue(QMStateType::Electron), "U_xN_xX_e");
  w(U_xN_xX_.getValue(QMStateType::Hole), "U_xN_xX_h");
  w(U_xN_xX_.getValue(QMStateType::Singlet), "U_xN_xX_s");
  w(U_xN_xX_.getValue(QMStateType::Triplet), "U_xN_xX_t");
  w(site_eng_.getValue(QMStateType::Electron), "site_eng_e");
  w(site_eng_.getValue(QMStateType::Hole), "site_eng_h");
  w(site_eng_.getValue(QMStateType::Singlet), "site_eng_s");
  w(site_eng_.getValue(QMStateType::Triplet), "site_eng_t");
}

void Segment::ReadFromCpt(CheckpointReader& r) {
  AtomContainer<Atom>::ReadFromCpt(r);
  r(molecule_ids_, "mol_ids");
  double value;
  r(value, "U_xX_nN_e");
  U_xX_nN_.setValue(value, QMStateType::Electron);
  r(value, "U_xX_nN_h");
  U_xX_nN_.setValue(value, QMStateType::Hole);
  r(value, "U_xX_nN_s");
  U_xX_nN_.setValue(value, QMStateType::Singlet);
  r(value, "U_xX_nN_t");
  U_xX_nN_.setValue(value, QMStateType::Triplet);
  r(value, "U_nX_nN_e");
  U_nX_nN_.setValue(value, QMStateType::Electron);
  r(value, "U_nX_nN_h");
  U_nX_nN_.setValue(value, QMStateType::Hole);
  r(value, "U_nX_nN_s");
  U_nX_nN_.setValue(value, QMStateType::Singlet);
  r(value, "U_nX_nN_t");
  U_nX_nN_.setValue(value, QMStateType::Triplet);
  r(value, "U_xN_xX_e");
  U_xN_xX_.setValue(value, QMStateType::Electron);
  r(value, "U_xN_xX_h");
  U_xN_xX_.setValue(value, QMStateType::Hole);
  r(value, "U_xN_xX_s");
  U_xN_xX_.setValue(value, QMStateType::Singlet);
  r(value, "U_xN_xX_t");
  U_xN_xX_.setValue(value, QMStateType::Triplet);
  r(value, "site_eng_e");
  site_eng_.setValue(value, QMStateType::Electron);
  r(value, "site_eng_h");
  site_eng_.setValue(value, QMStateType::Hole);
  r(value, "site_eng_s");
  site_eng_.setValue(value, QMStateType::Singlet);
  r(value, "site_eng_t");
  site_eng_.setValue(value, QMStateType::Triplet);
}

}  // namespace xtp
}  // namespace votca
