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

namespace votca {
namespace xtp {

QMPair::QMPair(int id, Segment *seg1, Segment *seg2,
               const Eigen::Vector3d &delta_R)
    : _id(id), _R(delta_R) {
  _segments.first = seg1;
  _segments.second = seg2;
}

Segment *QMPair::Seg2PbCopy() {
  Eigen::Vector3d r1 = _segments.first->getPos().toEigen();
  Eigen::Vector3d r2 = _segments.second->getPos().toEigen();

  // Check whether pair formed across periodic boundary
  if ((r2 - r1 - _R).norm() > 1e-8) {
    _ghost = std::unique_ptr<Segment>(new Segment(seg2));
    _ghost->TranslateBy(r1 - r2 + _R);
  }

  if (_ghost != nullptr) {
    return _ghost.get();
  } else {
    return _segments.second;
  }
}

}  // namespace xtp
}  // namespace votca