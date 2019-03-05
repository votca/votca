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

#include <votca/xtp/qmpair.h>
#include <votca/xtp/segment.h>

#include "votca/xtp/atomcontainer.h"

namespace votca {
namespace xtp {

QMPair::QMPair(int id, const Segment* seg1, const Segment* seg2,
               const Eigen::Vector3d& delta_R)
    : _id(id), _R(delta_R) {
  _segments.first = seg1;
  _segments.second = seg2;
}

const Segment* QMPair::Seg2PbCopy() {
  const Eigen::Vector3d& r1 = _segments.first->getPos();
  const Eigen::Vector3d& r2 = _segments.second->getPos();

  // Check whether pair formed across periodic boundary
  if (_ghost == nullptr && (r2 - r1 - _R).norm() > 1e-8) {
    _ghost = std::unique_ptr<Segment>(new Segment(*(_segments.second)));
    _ghost->Translate(r1 - r2 + _R);
  }

  if (_ghost != nullptr) {
    return _ghost.get();
  } else {
    return _segments.second;
  }
}

void QMPair::WriteToCpt(CheckpointWriter& w) const {
  w(_id, "index");
  w(_R, "delta_R");
  w(get_name(_pair_type), "pairtype");
  _lambda0.WriteToCpt(w, "lambdaO");
  _Jeff2.WriteToCpt(w, "Jeff2");
  w(_segments.first->getId(), "Seg1Id");
  w(_segments.second->getId(), "Seg2Id");
}

void QMPair::ReadFromCpt(CheckpointReader& r,
                         const std::vector<Segment>& segments) {
  r(_id, "index");
  r(_R, "delta_R");
  _lambda0.ReadFromCpt(r, "lambdaO");
  _Jeff2.ReadFromCpt(r, "Jeff2");
  std::string type_enum;
  r(type_enum, "pairtype");
  _pair_type = QMPair::get_Enum(type_enum);
  int id1;
  r(id1, "Seg1Id");
  int id2;
  r(id2, "Seg2Id");
  _segments.first = &segments[id1];
  _segments.second = &segments[id2];
}

}  // namespace xtp
}  // namespace votca