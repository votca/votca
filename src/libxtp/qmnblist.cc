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

#include <iostream>
#include <votca/xtp/qmnblist.h>

#include "votca/xtp/checkpointwriter.h"

using namespace std;

namespace votca {
namespace xtp {

QMPair& QMNBList::Add(Segment* seg1, Segment* seg2, const Eigen::Vector3d& r) {
  assert(this->FindPair(seg1, seg2) != NULL &&
         "Critical bug: pair already exists");
  int id = this->size();
  QMPair* pair = new QMPair(id, seg1, seg2, r);
  this->AddPair(pair);
  return *pair;
}

void QMNBList::WriteToCpt(CheckpointWriter& w) const {
  w(size(), "size");
  for (const QMPair* pair : _pairs) {
    CheckpointWriter u = w.openChild("pair" + std::to_string(pair->getId()));
    pair->WriteToCpt(u);
  }
}

void QMNBList::ReadFromCpt(CheckpointReader& r,
                           const std::vector<Segment>& segments) {
  Cleanup();
  int size = 0;
  r(size, "size");
  int count = r.getNumDataSets();
  if (count != size) {
    throw std::runtime_error(
        "QMNBList: size of neighborlist and number of saved pairs does not "
        "agree.");
  }
  for (int i = 0; i < size; i++) {
    CheckpointReader c = r.openChild("pair" + std::to_string(i));
    this->AddPair(new QMPair(c, segments));
  }
}

}  // namespace xtp
}  // namespace votca
