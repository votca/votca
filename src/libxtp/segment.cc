/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <votca/xtp/segment.h>

#include "votca/xtp/atomcontainer.h"

namespace votca {
namespace xtp {

double Segment::getApproxSize() const {
  std::pair<Eigen::Vector3d, Eigen::Vector3d> minmax = CalcSpatialMinMax();
  return (minmax.first - minmax.second).norm();
}

const Atom* Segment::getAtom(int id) const {

  for (const Atom& atom : *this) {
    if (atom.getId() == id) {
      return &atom;
    }
  }
  return nullptr;
}

void Segment::WriteToCpt(CheckpointWriter& w) const {
  AtomContainer<Atom>::WriteToCpt(w);
  w(_molecule_ids, "mol_ids");

  w(_U_xX_nN.getValue(QMStateType::Electron), "U_xX_nN_e");
  w(_U_xX_nN.getValue(QMStateType::Hole), "U_xX_nN_h");
  w(_U_xX_nN.getValue(QMStateType::Singlet), "U_xX_nN_s");
  w(_U_xX_nN.getValue(QMStateType::Triplet), "U_xX_nN_t");
  w(_U_nX_nN.getValue(QMStateType::Electron), "U_nX_nN_e");
  w(_U_nX_nN.getValue(QMStateType::Hole), "U_nX_nN_h");
  w(_U_nX_nN.getValue(QMStateType::Singlet), "U_nX_nN_s");
  w(_U_nX_nN.getValue(QMStateType::Triplet), "U_nX_nN_t");
  w(_U_xN_xX.getValue(QMStateType::Electron), "U_xN_xX_e");
  w(_U_xN_xX.getValue(QMStateType::Hole), "U_xN_xX_h");
  w(_U_xN_xX.getValue(QMStateType::Singlet), "U_xN_xX_s");
  w(_U_xN_xX.getValue(QMStateType::Triplet), "U_xN_xX_t");
  w(_site_eng.getValue(QMStateType::Electron), "site_eng_e");
  w(_site_eng.getValue(QMStateType::Hole), "site_eng_h");
  w(_site_eng.getValue(QMStateType::Singlet), "site_eng_s");
  w(_site_eng.getValue(QMStateType::Triplet), "site_eng_t");
}

void Segment::ReadFromCpt(CheckpointReader& r) {
  AtomContainer<Atom>::ReadFromCpt(r);
  r(_molecule_ids, "mol_ids");
  double value;
  r(value, "U_xX_nN_e");
  _U_xX_nN.setValue(value, QMStateType::Electron);
  r(value, "U_xX_nN_h");
  _U_xX_nN.setValue(value, QMStateType::Hole);
  r(value, "U_xX_nN_s");
  _U_xX_nN.setValue(value, QMStateType::Singlet);
  r(value, "U_xX_nN_t");
  _U_xX_nN.setValue(value, QMStateType::Triplet);
  r(value, "U_nX_nN_e");
  _U_nX_nN.setValue(value, QMStateType::Electron);
  r(value, "U_nX_nN_h");
  _U_nX_nN.setValue(value, QMStateType::Hole);
  r(value, "U_nX_nN_s");
  _U_nX_nN.setValue(value, QMStateType::Singlet);
  r(value, "U_nX_nN_t");
  _U_nX_nN.setValue(value, QMStateType::Triplet);
  r(value, "U_xN_xX_e");
  _U_xN_xX.setValue(value, QMStateType::Electron);
  r(value, "U_xN_xX_h");
  _U_xN_xX.setValue(value, QMStateType::Hole);
  r(value, "U_xN_xX_s");
  _U_xN_xX.setValue(value, QMStateType::Singlet);
  r(value, "U_xN_xX_t");
  _U_xN_xX.setValue(value, QMStateType::Triplet);
  r(value, "site_eng_e");
  _site_eng.setValue(value, QMStateType::Electron);
  r(value, "site_eng_h");
  _site_eng.setValue(value, QMStateType::Hole);
  r(value, "site_eng_s");
  _site_eng.setValue(value, QMStateType::Singlet);
  r(value, "site_eng_t");
  _site_eng.setValue(value, QMStateType::Triplet);
}

}  // namespace xtp
}  // namespace votca
