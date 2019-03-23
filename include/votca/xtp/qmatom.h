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

#ifndef VOTCA_XTP_QMATOM_H
#define VOTCA_XTP_QMATOM_H

#include <votca/tools/elements.h>
#include <votca/xtp/checkpoint.h>

namespace votca {
namespace xtp {

/**
 *    \brief container for QM atoms
 *
 *    Stores atom type, coordinates, charge
 */
class QMAtom {
  friend class AOBasis;

 public:
  QMAtom(int index, std::string element, Eigen::Vector3d pos)
      : _index(index),
        _element(element),
        _pos(pos),
        _nuccharge(0),
        _ecpcharge(0) {
    tools::Elements elements;
    _nuccharge = elements.getNucCrg(_element);
  }

  QMAtom(const CheckpointReader& r) { ReadFromCpt(r); }
  QMAtom(CptTable& table, const std::size_t& idx) { ReadFromCpt(table, idx); }

  const Eigen::Vector3d& getPos() const { return _pos; }

  void Translate(const Eigen::Vector3d& shift) { _pos += shift; }

  void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& refPos) {
    Eigen::Vector3d dir=_pos - refPos;
    dir = R * dir;
    _pos = refPos + dir;  // Rotated Position
  }

  void setPos(const Eigen::Vector3d& position) { _pos = position; }

  const std::string& getElement() const { return _element; }

  int getId() const { return _index; }

  int getNuccharge() const { return _nuccharge - _ecpcharge; }

  std::string identify() const { return "qmatom"; }

  friend std::ostream &operator<<(std::ostream &out, const QMAtom& atom) {
    out <<atom.getId()<<" "<<atom.getElement();
    out <<" "<<atom.getPos().x()<<","
	    <<atom.getPos().y()<<","<<atom.getPos().z()<<" "<<atom.getNuccharge()<<"\n";
    return out;
  }

 private:
  int _index;
  std::string _element;
  Eigen::Vector3d _pos;  // Bohr
  int _nuccharge;        // nuc charge is set in aobasis fill and ecpfill
  int _ecpcharge;

 public:
  struct data{
      int index;
      std::string element;
      double x;
      double y;
      double z;
      int nuccharge;
      int ecpcharge;
  };

  void SetupCptTable(CptTable& table) const {
      table.addCol(_index, "index", HOFFSET(data, index));
      table.addCol(_element, "element", HOFFSET(data, element));
      table.addCol(_pos[0], "pos.x", HOFFSET(data, x));
      table.addCol(_pos[1],"pos.y", HOFFSET(data, y));
      table.addCol(_pos[2],"pos.z", HOFFSET(data, z));
      table.addCol(_nuccharge, "nuccharge", HOFFSET(data, nuccharge));
      table.addCol(_ecpcharge, "ecpcharge", HOFFSET(data, ecpcharge));
  }

  void WriteToCpt(CptTable& table, const std::size_t& idx) const{
      data d;
      d.index     = _index;
      d.element   = _element;
      d.x         = _pos[0];
      d.y         = _pos[1];
      d.z         = _pos[2];
      d.nuccharge = _nuccharge;
      d.ecpcharge = _ecpcharge;

      table.writeToRow(&d, idx);
  };


  void WriteToCpt(const CheckpointWriter& w) const {
    w(_index, "index");
    w(_element, "element");
    w(_pos, "pos");
    w(_nuccharge, "nuccharge");
    w(_ecpcharge, "ecpcharge");
  }

  void ReadFromCpt(CptTable& table, const std::size_t& idx){
      data d;
      d.element = std::string("something really very long");
      table.readFromRow(&d, idx);

      _index     = d.index;
      _element   = std::string(d.element.c_str());
      _pos[0]    = d.x;
      _pos[1]    = d.y;
      _pos[2]    = d.z;
      _nuccharge = d.nuccharge;
      _ecpcharge = d.ecpcharge;
  }

  void ReadFromCpt(const CheckpointReader& r) {
    r(_index, "index");
    r(_element, "element");
    r(_pos, "pos");
    r(_nuccharge, "nuccharge");
    r(_ecpcharge, "ecpcharge");
  }
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMATOM_H
