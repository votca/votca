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

// Local VOTCA includes
#include "votca/xtp/qmpair.h"
#include "votca/xtp/atomcontainer.h"
#include "votca/xtp/segment.h"

namespace votca {
namespace xtp {

QMPair::QMPair(Index id, const Segment* seg1, const Segment* seg2,
               const Eigen::Vector3d& delta_R)
    : _id(id), _R(delta_R) {
  _segments.first = seg1;
  _segments.second = seg2;
}

Segment QMPair::Seg2PbCopy() const {
  const Eigen::Vector3d& r1 = _segments.first->getPos();
  const Eigen::Vector3d& r2 = _segments.second->getPos();
  Segment seg2pbc = Segment(*(_segments.second));
  // Check whether pair formed across periodic boundary
  if ((r2 - r1 - _R).norm() > 1e-8) {
    seg2pbc.Translate(r1 - r2 + _R);
  }

  return seg2pbc;
}

void QMPair::SetupCptTable(CptTable& table) {
  table.addCol<Index>("index", HOFFSET(data, id));
  table.addCol<Index>("Seg1Id", HOFFSET(data, Seg1Id));
  table.addCol<Index>("Seg2Id", HOFFSET(data, Seg2Id));

  table.addCol<double>("delta_Rx", HOFFSET(data, RX));
  table.addCol<double>("delta_Ry", HOFFSET(data, RY));
  table.addCol<double>("delta_Rz", HOFFSET(data, RZ));
  table.addCol<std::string>("pair_type", HOFFSET(data, pair_type));

  table.addCol<double>("lambda0e", HOFFSET(data, lambda0e));
  table.addCol<double>("lambda0h", HOFFSET(data, lambda0h));
  table.addCol<double>("lambda0s", HOFFSET(data, lambda0s));
  table.addCol<double>("lambda0t", HOFFSET(data, lambda0t));

  table.addCol<double>("jeff2e", HOFFSET(data, jeff2e));
  table.addCol<double>("jeff2h", HOFFSET(data, jeff2h));
  table.addCol<double>("jeff2s", HOFFSET(data, jeff2s));
  table.addCol<double>("jeff2t", HOFFSET(data, jeff2t));
}

void QMPair::WriteData(data& d) const {
  d.id = _id;
  d.Seg1Id = _segments.first->getId();
  d.Seg2Id = _segments.second->getId();
  d.RX = _R[0];
  d.RY = _R[1];
  d.RZ = _R[2];
  std::string ptype = get_name(_pair_type);
  d.pair_type = new char[ptype.length() + 1];
  strcpy(d.pair_type, ptype.c_str());

  d.lambda0e = _lambda0.getValue(QMStateType::Electron);
  d.lambda0h = _lambda0.getValue(QMStateType::Hole);
  d.lambda0s = _lambda0.getValue(QMStateType::Singlet);
  d.lambda0t = _lambda0.getValue(QMStateType::Triplet);

  d.jeff2e = _Jeff2.getValue(QMStateType::Electron);
  d.jeff2h = _Jeff2.getValue(QMStateType::Hole);
  d.jeff2s = _Jeff2.getValue(QMStateType::Singlet);
  d.jeff2t = _Jeff2.getValue(QMStateType::Triplet);
}

void QMPair::ReadData(const data& d, const std::vector<Segment>& segments) {
  _id = d.id;
  _R[0] = d.RX;
  _R[1] = d.RY;
  _R[2] = d.RZ;

  std::string type_enum = std::string(d.pair_type);
  _pair_type = QMPair::get_Enum(type_enum);
  free(d.pair_type);

  _lambda0.setValue(d.lambda0e, QMStateType::Electron);
  _lambda0.setValue(d.lambda0h, QMStateType::Hole);
  _lambda0.setValue(d.lambda0s, QMStateType::Singlet);
  _lambda0.setValue(d.lambda0t, QMStateType::Triplet);

  _Jeff2.setValue(d.jeff2e, QMStateType::Electron);
  _Jeff2.setValue(d.jeff2h, QMStateType::Hole);
  _Jeff2.setValue(d.jeff2s, QMStateType::Singlet);
  _Jeff2.setValue(d.jeff2t, QMStateType::Triplet);

  _segments.first = &segments[d.Seg1Id];
  _segments.second = &segments[d.Seg2Id];
}

}  // namespace xtp
}  // namespace votca
