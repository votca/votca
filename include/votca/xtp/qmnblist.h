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

#ifndef VOTCA_XTP_QMNBList_H
#define VOTCA_XTP_QMNBList_H

#include <votca/csg/pairlist.h>
#include <votca/xtp/qmpair.h>

namespace votca {
namespace xtp {

class QMNBList : public csg::PairList<Segment*, QMPair> {
 public:
  QMNBList(){};
  ~QMNBList() { csg::PairList<Segment*, QMPair>::Cleanup(); }

  QMPair& Add(Segment* seg1, Segment* seg2, const Eigen::Vector3d& r);

  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r, const std::vector<Segment>& segments);

 protected:
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMNBLIST_H
