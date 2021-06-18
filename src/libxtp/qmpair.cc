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
    : id_(id), R_(delta_R) {
  segments_.first = seg1;
  segments_.second = seg2;
}

Segment QMPair::Seg2PbCopy() const {
  const Eigen::Vector3d& r1 = segments_.first->getPos();
  const Eigen::Vector3d& r2 = segments_.second->getPos();
  Segment seg2pbc = Segment(*(segments_.second));
  // Check whether pair formed across periodic boundary
  if ((r2 - r1 - R_).norm() > 1e-8) {
    seg2pbc.Translate(r1 - r2 + R_);
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
  d.id = id_;
  d.Seg1Id = segments_.first->getId();
  d.Seg2Id = segments_.second->getId();
  d.RX = R_[0];
  d.RY = R_[1];
  d.RZ = R_[2];
  std::string ptype = get_name(pair_type_);
  d.pair_type = new char[ptype.length() + 1];
  strcpy(d.pair_type, ptype.c_str());

  d.lambda0e = lambda0_.getValue(QMStateType::Electron);
  d.lambda0h = lambda0_.getValue(QMStateType::Hole);
  d.lambda0s = lambda0_.getValue(QMStateType::Singlet);
  d.lambda0t = lambda0_.getValue(QMStateType::Triplet);

  d.jeff2e = Jeff2_.getValue(QMStateType::Electron);
  d.jeff2h = Jeff2_.getValue(QMStateType::Hole);
  d.jeff2s = Jeff2_.getValue(QMStateType::Singlet);
  d.jeff2t = Jeff2_.getValue(QMStateType::Triplet);
}

void QMPair::ReadData(const data& d, const std::vector<Segment>& segments) {
  id_ = d.id;
  R_[0] = d.RX;
  R_[1] = d.RY;
  R_[2] = d.RZ;

  std::string type_enum = std::string(d.pair_type);
  pair_type_ = QMPair::get_Enum(type_enum);
  free(d.pair_type);

  lambda0_.setValue(d.lambda0e, QMStateType::Electron);
  lambda0_.setValue(d.lambda0h, QMStateType::Hole);
  lambda0_.setValue(d.lambda0s, QMStateType::Singlet);
  lambda0_.setValue(d.lambda0t, QMStateType::Triplet);

  Jeff2_.setValue(d.jeff2e, QMStateType::Electron);
  Jeff2_.setValue(d.jeff2h, QMStateType::Hole);
  Jeff2_.setValue(d.jeff2s, QMStateType::Singlet);
  Jeff2_.setValue(d.jeff2t, QMStateType::Triplet);

  segments_.first = &segments[d.Seg1Id];
  segments_.second = &segments[d.Seg2Id];
}

}  // namespace xtp
}  // namespace votca
