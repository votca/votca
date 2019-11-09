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
#include "votca/xtp/atom.h"
#include <votca/tools/elements.h>

namespace votca {
namespace xtp {

Atom::Atom(Index resnr, std::string md_atom_name, Index atom_id,
           Eigen::Vector3d pos, std::string type)
    : _id(atom_id), _name(md_atom_name), _resnr(resnr), _pos(pos) {

  std::string elename = GetElementFromString(md_atom_name);
  std::string eletype = GetElementFromString(type);
  tools::Elements ele;
  bool found_element_name = true;
  bool found_element_type = true;
  try {
    ele.getMass(elename);
  } catch (std::runtime_error&) {
    found_element_name = false;
  }

  try {
    ele.getMass(eletype);
  } catch (std::runtime_error&) {
    found_element_type = false;
  }

  if (found_element_name && found_element_type) {
    if (elename != eletype) {
      throw std::runtime_error("Elements " + elename + " and" + eletype +
                               " from atom name: " + md_atom_name +
                               " and atom type:" + type + " do not match.");
    }
    _element = elename;
  } else if (found_element_name) {
    _element = elename;
  } else if (found_element_type) {
    _element = elename;
  } else {
    throw std::runtime_error("Could not get Element from atom name:" +
                             md_atom_name + " or atom type:" + type);
  }
}

Atom::Atom(Index atom_id, std::string element, Eigen::Vector3d pos)
    : Atom(-1, element, atom_id, pos, element) {}

std::string Atom::GetElementFromString(const std::string& MDName) {
  std::string element = MDName.substr(0, 1);

  if (MDName.size() > 1) {
    if (std::islower(MDName[1])) {
      element += MDName[1];
    }
  }
  return element;
}

void Atom::Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& refPos) {
  Eigen::Vector3d dir = _pos - refPos;
  dir = R * dir;
  _pos = refPos + dir;  // Rotated Position
}

void Atom::SetupCptTable(CptTable& table) const {
  table.addCol(_id, "index", HOFFSET(data, id));
  table.addCol(_element, "element", HOFFSET(data, element));
  table.addCol(_name, "name", HOFFSET(data, name));
  table.addCol(_pos[0], "pos.x", HOFFSET(data, x));
  table.addCol(_pos[1], "pos.y", HOFFSET(data, y));
  table.addCol(_pos[2], "pos.z", HOFFSET(data, z));
  table.addCol(_resnr, "resnr", HOFFSET(data, resnr));
}

void Atom::WriteData(data& d) const {
  d.id = _id;
  d.element = const_cast<char*>(_element.c_str());
  d.name = const_cast<char*>(_name.c_str());
  d.x = _pos[0];
  d.y = _pos[1];
  d.z = _pos[2];
  d.resnr = _resnr;
}

void Atom::ReadData(const data& d) {
  _id = d.id;
  _element = std::string(d.element);
  free(d.element);
  _name = std::string(d.name);
  free(d.name);
  _pos[0] = d.x;
  _pos[2] = d.z;
  _pos[1] = d.y;
  _resnr = d.resnr;
}
}  // namespace xtp
}  // namespace votca
