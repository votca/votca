/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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

#ifndef __VOTCA_XTP_ATOM_H
#define	__VOTCA_XTP_ATOM_H

#include <string>
#include <assert.h>
#include <votca/tools/types.h>
#include <votca/tools/vec.h>
#include <votca/tools/property.h>
#include <votca/tools/matrix.h>




namespace votca { namespace xtp {
using namespace votca::tools;

using namespace std;

class Topology;
class Molecule;
class Segment;
class Fragment;

/**
    \brief information about an atom
 
    The Atom class stores atom id, name, type, mass, charge, residue number
    
*/
class Atom 
{
public:   

    Atom(Molecule *owner,
         string residue_name,   int resnr,
         string md_atom_name,   int md_atom_id,
         bool hasQMPart,        int qm_atom_id,
         vec qmPos,             string element,
         double weight)
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
	 _element(element)  { }
    
    Atom(int atom_id,   string atom_name)
       : _id(atom_id),  _name(atom_name),
         _hasQM(false), _qmId(-1) { }

    Atom(Atom *stencil)
       : _id(stencil->getId()),
         _name(stencil->getName()+"_ghost"),
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

    Atom() { };
   ~Atom() { _Q.clear(); }

    const int     &getId() const { return _id; }
    const string  &getName() const { return _name; }
    const string  &getType() const { return _type; }
    const int     &getResnr() const { return _resnr; }

    inline void setTopology(Topology *container) { _top = container; }
    inline void setMolecule(Molecule *container) { _mol = container; }
    inline void setSegment(Segment *container)   { _seg = container; }
    inline void setFragment(Fragment *container) { _frag = container; }

    Topology *getTopology() { return _top; }
    Molecule *getMolecule() { return _mol; }
    Segment  *getSegment() { return _seg; }
    Fragment *getFragment() { return _frag; }

    inline void setResnr(const int &resnr) { _resnr = resnr; }
    inline void setResname(const string &resname) { _resname = resname; }
    inline void setWeight(const double &weight) { _weight = weight; }
    inline void setQMPart(const int &qmid, vec qmPos);
    inline void setQMPos(const vec &qmPos) { _qmPos = qmPos; }
    inline void setElement(const string &element) { _element = element; }
    inline void TranslateBy(const vec &shift) { _pos = _pos + shift;  }
    
    inline const int    &getResnr() { return _resnr; }
    inline const string &getResname() { return _resname; }
    inline const double &getWeight() { return _weight; }
    inline const int    &getQMId() { return _qmId; }
    inline const vec    &getQMPos() { return _qmPos; }
    inline const string &getElement() { return _element; }

    inline const double &getQ(int state) { return _Q.at(state); }
    inline const double &getQ() { return _q->second; }
    inline void          setQ( map<int, double> Q) { _Q = Q; }
    void                 chrg(int state) { _q = _Q.find(state); }

    inline void          setPTensor(matrix &ptensor) { _ptensor = ptensor; }
    const matrix        &getPTensor() { return _ptensor; }

    /**
     * get the position of the atom
     * \return atom position
     */
    const vec &getPos() const;
    /**
     * set the position of the atom
     * \param r atom position
     */
    void setPos(const vec &r);
    /**
     * direct access (read/write) to the position of the atom
     * \return reference to position 
     */
    vec &Pos() { return _pos; }
    /** does this configuration store positions? */
    bool HasPos() {return _bPos; }           
    /** dose the bead store a position */
    void HasPos(bool b);


    bool HasQMPart() { return _hasQM; }
    /**
     * molecule the bead belongs to
     * \return Molecule object
     */







protected:
    int         _id;
    string      _name;

    Topology   *_top;
    Molecule   *_mol;
    Segment    *_seg;
    Fragment   *_frag;    

    string      _type;
    int         _resnr;
    string      _resname;
    double      _weight;
    vec         _pos;
    bool        _bPos;

    bool        _hasQM;
    int         _qmId;
    vec         _qmPos;
    string      _element;

    // charge state of segment => partial charge
    map<int, double> _Q;    
    map<int, double> ::iterator _q;
    matrix _ptensor;
        
};

inline void Atom::setPos(const vec &r) {
    _bPos=true;
    _pos = r;
}

inline const vec &Atom::getPos() const {
    assert(_bPos);
    return _pos;
}

inline void Atom::HasPos(bool b) {
    _bPos=b;
}

inline void Atom::setQMPart(const int &qmid, vec qmPos) {
    if (qmid > -1) {
        _hasQM = true;
        _qmId = qmid;
        _qmPos = qmPos;
    }
    else {
        _hasQM = false;
        _qmId = -1;
    }
}

}}

#endif	/* __VOTCA_XTP_ATOM_H */

