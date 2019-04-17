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
#include "votca/xtp/qmatom.h"
#include <iostream>

namespace votca {
namespace xtp {

QMAtom::QMAtom(int index, std::string element, Eigen::Vector3d pos)
    : _index(index), _element(element), _pos(pos) {
  tools::Elements elements;
  _nuccharge = elements.getNucCrg(_element);
}

QMAtom::QMAtom(data& d) { ReadData(d); }

void QMAtom::Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& refPos) {
  Eigen::Vector3d dir = _pos - refPos;
  dir = R * dir;
  _pos = refPos + dir;  // Rotated Position
}

void QMAtom::SetupCptTable(CptTable& table) const {
  table.addCol(_index, "index", HOFFSET(data, index));
  table.addCol(_element, "element", HOFFSET(data, element));
  table.addCol(_pos[0], "pos.x", HOFFSET(data, x));
  table.addCol(_pos[1], "pos.y", HOFFSET(data, y));
  table.addCol(_pos[2], "pos.z", HOFFSET(data, z));
  table.addCol(_nuccharge, "nuccharge", HOFFSET(data, nuccharge));
  table.addCol(_ecpcharge, "ecpcharge", HOFFSET(data, ecpcharge));
}

void QMAtom::WriteToCpt(CptTable& table, const std::size_t& idx) const {
  data d;
  d.index = _index;
  std::strcpy(d.element, _element.c_str());
  d.x = _pos[0];
  d.y = _pos[1];
  d.z = _pos[2];
  d.nuccharge = _nuccharge;
  d.ecpcharge = _ecpcharge;

  table.writeToRow(&d, idx);
};

void QMAtom::WriteData(data& d) const {
  d.index = _index;
  std::strcpy(d.element, _element.c_str());
  d.x = _pos[0];
  d.y = _pos[1];
  d.z = _pos[2];
  d.nuccharge = _nuccharge;
  d.ecpcharge = _ecpcharge;
}

void QMAtom::ReadFromCpt(CptTable& table, const std::size_t& idx) {
  data d;
  table.readFromRow(&d, idx);
  _element = std::string(d.element);
  // free(d.element);
  _index = d.index;
  _pos[0] = d.x;
  _pos[1] = d.y;
  _pos[2] = d.z;
  _nuccharge = d.nuccharge;
  _ecpcharge = d.ecpcharge;
}

void QMAtom::ReadData(data& d) {
  _element = std::string(d.element);
  _index = d.index;
  _pos[0] = d.x;
  _pos[1] = d.y;
  _pos[2] = d.z;
  _nuccharge = d.nuccharge;
  _ecpcharge = d.ecpcharge;
}

}  // namespace xtp
}  // namespace votca
