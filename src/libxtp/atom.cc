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

// VOTCA includes
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/atom.h"
#include "votca/xtp/checkpointtable.h"

namespace votca {
namespace xtp {

Atom::Atom(Index resnr, std::string md_atom_name, Index atom_id,
           Eigen::Vector3d pos, std::string type)
    : id_(atom_id), name_(md_atom_name), resnr_(resnr), pos_(pos) {

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
    element_ = elename;
  } else if (found_element_name) {
    element_ = elename;
  } else if (found_element_type) {
    element_ = elename;
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
  Eigen::Vector3d dir = pos_ - refPos;
  dir = R * dir;
  pos_ = refPos + dir;  // Rotated Position
}

void Atom::SetupCptTable(CptTable& table) {
  table.addCol<Index>("index", HOFFSET(data, id));
  table.addCol<std::string>("element", HOFFSET(data, element));
  table.addCol<std::string>("name", HOFFSET(data, name));
  table.addCol<double>("pos.x", HOFFSET(data, x));
  table.addCol<double>("pos.y", HOFFSET(data, y));
  table.addCol<double>("pos.z", HOFFSET(data, z));
  table.addCol<Index>("resnr", HOFFSET(data, resnr));
}

void Atom::WriteData(data& d) const {
  d.id = id_;
  d.element = const_cast<char*>(element_.c_str());
  d.name = const_cast<char*>(name_.c_str());
  d.x = pos_[0];
  d.y = pos_[1];
  d.z = pos_[2];
  d.resnr = resnr_;
}

void Atom::ReadData(const data& d) {
  id_ = d.id;
  element_ = std::string(d.element);
  free(d.element);
  name_ = std::string(d.name);
  free(d.name);
  pos_[0] = d.x;
  pos_[2] = d.z;
  pos_[1] = d.y;
  resnr_ = d.resnr;
}
}  // namespace xtp
}  // namespace votca
