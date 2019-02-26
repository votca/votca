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
/// For earlier commit history see ctp commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef __VOTCA_XTP_ATOM_H
#define __VOTCA_XTP_ATOM_H

#include <exception>
#include <map>
#include <string>
#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>

namespace votca {
namespace xtp {

class Topology;
class Molecule;
class Segment;
class Fragment;

/**
    \brief information about an atom

    The Atom class stores atom id, name, type, mass, charge, residue number

*/
class Atom {
 public:
  Atom(Molecule *owner, std::string residue_name, int resnr,
       std::string md_atom_name, int md_atom_id, bool hasQMPart, int qm_atom_id,
       votca::tools::vec qmPos, std::string element, double weight)
      : _id(md_atom_id),
        _name(md_atom_name),
        _mol(owner),
        _resnr(resnr),
        _resname(residue_name),
        _weight(weight),
        _bPos(false),
        _hasQM(hasQMPart),
        _qmId(qm_atom_id),
        _qmPos(qmPos),
        _element(element) {}

  Atom(int atom_id, std::string atom_name)
      : _id(atom_id), _name(atom_name), _hasQM(false), _qmId(-1) {}

  // TODO This should be replaced from a constructor to an overloaded = operator
  Atom(Atom *stencil)
      : _id(stencil->getId()),
        _name(stencil->getName() + "_ghost"),
        _top(NULL),
        _mol(NULL),
        _resnr(stencil->getResnr()),
        _resname(stencil->getResname()),
        _weight(stencil->getWeight()),
        _pos(stencil->getPos()),
        _bPos(true),
        _hasQM(stencil->HasQMPart()),
        _qmId(stencil->getQMId()),
        _qmPos(stencil->getQMPos()),
        _element(stencil->getElement()) {}

  Atom(){};
  ~Atom() { _Q.clear(); }

  const int &getId() const { return _id; }
  const std::string &getName() const { return _name; }
  const std::string &getType() const { return _type; }
  const int &getResnr() const { return _resnr; }

  inline void setTopology(Topology *container) { _top = container; }
  inline void setMolecule(Molecule *container) { _mol = container; }
  inline void setSegment(Segment *container) { _seg = container; }
  inline void setFragment(Fragment *container) { _frag = container; }

  Topology *getTopology() { return _top; }
  Molecule *getMolecule() { return _mol; }
  Segment *getSegment() { return _seg; }
  Fragment *getFragment() { return _frag; }

  inline void setResnr(const int &resnr) { _resnr = resnr; }
  inline void setResname(const std::string &resname) { _resname = resname; }
  inline void setWeight(const double &weight) { _weight = weight; }
  inline void setQMPart(const int &qmid, votca::tools::vec qmPos);
  inline void setQMPos(const votca::tools::vec &qmPos) { _qmPos = qmPos; }
  inline void setElement(const std::string &element) { _element = element; }
  inline void TranslateBy(const votca::tools::vec &shift) {
    _pos = _pos + shift;
  }

  inline const int &getResnr() { return _resnr; }
  inline const std::string &getResname() { return _resname; }
  inline const double &getWeight() { return _weight; }
  inline const int &getQMId() { return _qmId; }
  inline const votca::tools::vec &getQMPos() { return _qmPos; }
  inline const std::string &getElement() { return _element; }

  inline const double &getQ(int state) { return _Q.at(state); }
  inline const double &getQ() { return _q->second; }
  inline void setQ(std::map<int, double> Q) { _Q = Q; }
  void chrg(int state) { _q = _Q.find(state); }

  inline void setPTensor(votca::tools::matrix &ptensor) { _ptensor = ptensor; }
  const votca::tools::matrix &getPTensor() { return _ptensor; }

  /**
   * get the position of the atom
   * \return atom position
   */
  const votca::tools::vec &getPos() const;
  /**
   * set the position of the atom
   * \param r atom position
   */
  void setPos(const votca::tools::vec &r);
  /**
   * direct access (read/write) to the position of the atom
   * \return reference to position
   */
  votca::tools::vec &Pos() { return _pos; }
  /** does this configuration store positions? */
  bool HasPos() { return _bPos; }
  /** dose the bead store a position */
  void HasPos(bool b);

  bool HasQMPart() { return _hasQM; }
  /**
   * molecule the bead belongs to
   * \return Molecule object
   */

 protected:
  int _id;
  std::string _name;

  Topology *_top;
  Molecule *_mol;
  Segment *_seg;
  Fragment *_frag;

  std::string _type;
  int _resnr;
  std::string _resname;
  double _weight;
  votca::tools::vec _pos;
  bool _bPos;

  bool _hasQM;
  int _qmId;
  votca::tools::vec _qmPos;
  std::string _element;

  // charge state of segment => partial charge
  std::map<int, double> _Q;
  std::map<int, double>::iterator _q;
  votca::tools::matrix _ptensor;
};

inline void Atom::setPos(const votca::tools::vec &r) {
  _bPos = true;
  _pos = r;
}

inline const votca::tools::vec &Atom::getPos() const {
  if (!_bPos) throw std::runtime_error("Position has not yet been set");
  return _pos;
}

inline void Atom::HasPos(bool b) { _bPos = b; }

inline void Atom::setQMPart(const int &qmid, votca::tools::vec qmPos) {
  if (qmid > -1) {
    _hasQM = true;
    _qmId = qmid;
    _qmPos = qmPos;
  } else {
    _hasQM = false;
    _qmId = -1;
  }
}
}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_ATOM_H */
