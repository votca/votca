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
#ifndef VOTCA_XTP_QMNBLIST_H
#define VOTCA_XTP_QMNBLIST_H

// VOTCA includes
#include <votca/csg/pairlist.h>

// Local VOTCA includes
#include "qmpair.h"

namespace votca {
namespace xtp {

class QMNBList : public csg::PairList<const Segment*, QMPair> {
 public:
  QMNBList() = default;
  ~QMNBList() override { csg::PairList<const Segment*, QMPair>::Cleanup(); }

  QMPair& Add(const Segment& seg1, const Segment& seg2,
              const Eigen::Vector3d& r);

  template <class Compare>
  void sortAndReindex(Compare comp);

  const QMPair* operator[](Index index) const { return pairs_[index]; }
  QMPair* operator[](Index index) { return pairs_[index]; }

  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r, const std::vector<Segment>& segments);

 protected:
};

template <class Compare>
inline void QMNBList::sortAndReindex(Compare comp) {
  std::sort(pairs_.begin(), pairs_.end(), comp);

  for (Index i = 0; i < Index(pairs_.size()); i++) {
    pairs_[i]->setId(i);
  }
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMNBLIST_H
