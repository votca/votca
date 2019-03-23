/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

namespace votca {
namespace xtp {
    
     Atom::Atom(int resnr, std::string md_atom_name, int atom_id, Eigen::Vector3d pos)
      : _id(atom_id), _name(md_atom_name), _resnr(resnr), _pos(pos) {
    _element=GetElementFromMDName(md_atom_name);
  }

      
Atom::Atom(int atom_id,std::string md_atom_name, Eigen::Vector3d pos) : Atom(-1, md_atom_name, atom_id, pos)    {
      _element=GetElementFromMDName(md_atom_name);
    }

 std::string Atom::GetElementFromMDName(const std::string& MDName){
      std::string element = MDName.substr(0, 1);
        if (MDName.size() > 1){
            if (std::islower(MDName[1])){
                element += MDName[1];
            }
        }
      return element;
  }
 
 void Atom::Rotate(const Eigen::Matrix3d &R, const Eigen::Vector3d &refPos) {
   Eigen::Vector3d dir=_pos - refPos;
    dir = R * dir;
    _pos = refPos + dir;  // Rotated Position
  }
 
 
  void Atom::SetupCptTable(CptTable& table) const{
      table.addCol(_id, "index", HOFFSET(data, id));
      table.addCol(_element, "element", HOFFSET(data, element));
      table.addCol(_name, "name", HOFFSET(data, name));
      table.addCol(_pos[0], "pos.x", HOFFSET(data, x));
      table.addCol(_pos[1], "pos.y", HOFFSET(data, y));
      table.addCol(_pos[2], "pos.z", HOFFSET(data, z));
  }

  void Atom::WriteToCpt(CptTable& table, const std::size_t& idx) const{
      data d;
      d.id = _id;
      d.element = const_cast<char*>(_element.c_str());;
      d.name = const_cast<char*>(_name.c_str());;
      d.x = _pos[0];
      d.y = _pos[1];
      d.z = _pos[2];

      table.writeToRow(&d, idx);
  }

  void Atom::ReadFromCpt(CptTable table, const std::size_t& idx){
      data d;
      table.readFromRow(&d, idx);

      _id      = d.id;
      _element   = std::string(d.element);
      free(d.element);
      _name    = std::string(d.name);
       free(d.name);
      _pos[0]  = d.x;
      _pos[2]  = d.z;
      _pos[1]  = d.y;
  }
  

}  // namespace xtp
}  // namespace votca
