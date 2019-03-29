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

#ifndef VOTCA_XTP_ATOM_H
#define VOTCA_XTP_ATOM_H

#include <exception>
#include <map>
#include <string>
#include <votca/xtp/checkpointreader.h>
#include <votca/xtp/checkpointwriter.h>

namespace votca {
namespace xtp {
typedef std::pair<int,std::string> MD_atom_id;
/**
    \brief information about an atom

    The Atom class stores atom id, name, type,residue number

*/
class Atom {
 public:
    struct data{
        int id;
        char* element;
        char* name;
        double x;
        double y;
        double z;
    };
  Atom(int resnr, std::string md_atom_name, int atom_id, Eigen::Vector3d pos);

  Atom(int atom_id,std::string md_atom_name, Eigen::Vector3d pos);

  Atom(CptTable table, const std::size_t& idx ) { ReadFromCpt(table, idx); }

  Atom(data& d) { ReadData(d); }

  static std::string GetElementFromMDName(const std::string& MDName);

  int getId() const { return _id; }
  const std::string &getName() const { return _name; }
  std::string getElement() const { return _element; }

  int getResnr() const { return _resnr; }

  void setResnr(int resnr) { _resnr = resnr; }
  void Translate(const Eigen::Vector3d &shift) { _pos = _pos + shift; }

  void Rotate(const Eigen::Matrix3d &R, const Eigen::Vector3d &refPos);

  const Eigen::Vector3d &getPos() const { return _pos; }
  void setPos(const Eigen::Vector3d &r) { _pos = r; }

  std::string identify() const { return "atom"; }

  friend std::ostream &operator<<(std::ostream &out, const Atom& atom) {
    out <<atom.getId()<<" "<<atom.getName()<<" "<<atom.getElement()<<" "<<atom.getResnr();
    out <<" "<<atom.getPos().x()<<","<<
        atom.getPos().y()<<","<<atom.getPos().z()<<"\n";
    return out;
  }



  void SetupCptTable(CptTable& table) const;

  void WriteData(data& d) const;

  void WriteToCpt(CptTable& table, const std::size_t& idx) const;

  void ReadData(data& d);

  void ReadFromCpt(CptTable table, const std::size_t& idx);

 private:
  int _id = -1;
  std::string _name = "";

  std::string _element = "";
  int _resnr = -1;
  Eigen::Vector3d _pos = Eigen::Vector3d::Zero();
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ATOM_H
