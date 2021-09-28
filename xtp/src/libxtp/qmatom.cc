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
#include "votca/xtp/qmatom.h"
#include "votca/xtp/checkpointtable.h"
namespace votca {
namespace xtp {

QMAtom::QMAtom(Index index, std::string element, Eigen::Vector3d pos)
    : index_(index), element_(element), pos_(pos) {
  tools::Elements elements;
  nuccharge_ = elements.getNucCrg(element_);
}

QMAtom::QMAtom(const data& d) { ReadData(d); }

void QMAtom::Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& refPos) {
  Eigen::Vector3d dir = pos_ - refPos;
  dir = R * dir;
  pos_ = refPos + dir;  // Rotated Position
}

void QMAtom::SetupCptTable(CptTable& table) {
  table.addCol<Index>("index", HOFFSET(data, index));
  table.addCol<std::string>("element", HOFFSET(data, element));
  table.addCol<double>("posX", HOFFSET(data, x));
  table.addCol<double>("posY", HOFFSET(data, y));
  table.addCol<double>("posZ", HOFFSET(data, z));
  table.addCol<Index>("nuccharge", HOFFSET(data, nuccharge));
  table.addCol<Index>("ecpcharge", HOFFSET(data, ecpcharge));
}

void QMAtom::WriteData(data& d) const {
  d.index = index_;
  d.element = const_cast<char*>(element_.c_str());
  d.x = pos_[0];
  d.y = pos_[1];
  d.z = pos_[2];
  d.nuccharge = nuccharge_;
  d.ecpcharge = ecpcharge_;
}

void QMAtom::ReadData(const data& d) {
  element_ = std::string(d.element);
  free(d.element);
  index_ = d.index;
  pos_[0] = d.x;
  pos_[1] = d.y;
  pos_[2] = d.z;
  nuccharge_ = d.nuccharge;
  ecpcharge_ = d.ecpcharge;
}
}  // namespace xtp
}  // namespace votca
