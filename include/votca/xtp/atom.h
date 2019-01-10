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
/// For earlier commit history see ctp commit 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef VOTCA_XTP_ATOM_H
#define VOTCA_XTP_ATOM_H

#include <string>
#include <map>
#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>
#include <exception>
#include <map>

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
       tools::vec qmPos, std::string element, double weight)
      : _id(md_atom_id),
        _name(md_atom_name),
        _mol(owner),
        _resnr(resnr),
        _resname(residue_name),
        _weight(weight),
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
        _hasQM(stencil->HasQMPart()),
        _qmId(stencil->getQMId()),
        _qmPos(stencil->getQMPos()),
        _element(stencil->getElement()) {}

  Atom() {};

  const int &getId() const { return _id; }
  const std::string &getName() const { return _name; }
  const std::string &getType() const { return _type; }
  const int &getResnr() const { return _resnr; }

   void setTopology(Topology *container) { _top = container; }
   void setMolecule(Molecule *container) { _mol = container; }
   void setSegment(Segment *container) { _seg = container; }
   void setFragment(Fragment *container) { _frag = container; }

  Topology *getTopology() { return _top; }
  Molecule *getMolecule() { return _mol; }
  Segment *getSegment() { return _seg; }
  Fragment *getFragment() { return _frag; }

   void setResnr(int resnr) { _resnr = resnr; }
   void setResname(const std::string &resname) { _resname = resname; }
   void setWeight(double weight) { _weight = weight; }
   void setQMPart(int qmid, tools::vec qmPos);
   void setQMPos(const tools::vec &qmPos) { _qmPos = qmPos; }
   void setElement(const std::string &element) { _element = element; }
   void TranslateBy(const tools::vec &shift) {
    _pos = _pos + shift;
  }

   const int &getResnr() { return _resnr; }
   const std::string &getResname() { return _resname; }
   const double &getWeight() { return _weight; }
   const int &getQMId() { return _qmId; }
  const tools::vec &getQMPos() { return _qmPos; }
  const std::string &getElement() { return _element; }


  const tools::vec &getPos() const{return _pos;}
  void setPos(const tools::vec &r){_pos = r;}

  bool HasQMPart() { return _hasQM; }

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
  tools::vec _pos;

  bool _hasQM;
  int _qmId;
  tools::vec _qmPos;
  std::string _element;
};


inline void Atom::setQMPart(int qmid, tools::vec qmPos) {
  if (qmid > -1) {
    _hasQM = true;
    _qmId = qmid;
    _qmPos = qmPos;
  } else {
    _hasQM = false;
    _qmId = -1;
  }
}
}
}

#endif // VOTCA_XTP_ATOM_H 
