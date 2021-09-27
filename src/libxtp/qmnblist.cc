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
#include "votca/xtp/qmnblist.h"
#include "votca/xtp/checkpointwriter.h"

namespace votca {
namespace xtp {

QMPair& QMNBList::Add(const Segment& seg1, const Segment& seg2,
                      const Eigen::Vector3d& r) {
  assert(this->FindPair(&seg1, &seg2) == nullptr &&
         "Critical bug: pair already exists");
  Index id = this->size();
  QMPair* pair = new QMPair(id, &seg1, &seg2, r);
  this->AddPair(pair);
  return *pair;
}

void QMNBList::WriteToCpt(CheckpointWriter& w) const {
  Index size = this->size();
  w(size, "size");
  if (size == 0) {
    return;
  }
  std::vector<QMPair::data> dataVec(size);

  CptTable table = w.openTable<QMPair>("pairs", size);
  for (Index i = 0; i < size; i++) {
    (pairs_[i]->WriteData(dataVec[i]));
  }
  table.write(dataVec);
  for (QMPair::data data : dataVec) {
    delete[] data.pair_type;
  }
}

void QMNBList::ReadFromCpt(CheckpointReader& r,
                           const std::vector<Segment>& segments) {
  Cleanup();
  Index size = 0;
  r(size, "size");
  if (size == 0) {
    return;
  }
  CptTable table = r.openTable<QMPair>("pairs");
  std::vector<QMPair::data> dataVec(table.numRows());
  table.read(dataVec);

  for (const QMPair::data& data : dataVec) {
    this->AddPair(new QMPair(data, segments));
  }
}

}  // namespace xtp
}  // namespace votca
