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
/// For earlier commit history see ctp commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#include <votca/xtp/segment.h>

#include "votca/xtp/atomcontainer.h"

namespace votca {
namespace xtp {

double Segment::getApproxSize() const {
  if (!_has_approxsize || !this->PosIsValid()) {
    if (!PosIsValid()) {
      this->getPos();
    }
    std::pair<Eigen::Vector3d, Eigen::Vector3d> minmax = CalcSpatialMinMax();
    _approxsize = (minmax.first - minmax.second).norm();
    _has_approxsize = true;
  }
  return _approxsize;
}

void Segment::WriteToCpt(CheckpointWriter& w) const {
  AtomContainer<Atom>::WriteToCpt(w);
  w(_molecule_ids, "mol_ids");
  _U_xX_nN.WriteToCpt(w, "U_xX_nN");
  _U_nX_nN.WriteToCpt(w, "U_nX_nN");
  _U_xN_xX.WriteToCpt(w, "U_xN_xX");
  _eMpoles.WriteToCpt(w, "elStatic");
}

void Segment::ReadFromCpt(CheckpointReader& r) {
  AtomContainer<Atom>::ReadFromCpt(r);
  r(_molecule_ids, "mol_ids");
  _U_xX_nN.ReadFromCpt(r, "U_xX_nN");
  _U_nX_nN.ReadFromCpt(r, "U_nX_nN");
  _U_xN_xX.ReadFromCpt(r, "U_xN_xX");
  _eMpoles.ReadFromCpt(r, "elStatic");
}

}  // namespace xtp
}  // namespace votca
