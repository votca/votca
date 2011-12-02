/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __VOTCA_CTP_ATOM_H
#define	__VOTCA_CTP_ATOM_H

#include <string>
#include <assert.h>
#include <votca/tools/types.h>
#include <votca/tools/vec.h>
#include <votca/tools/property.h>


namespace votca { namespace ctp {
using namespace votca::tools;

using namespace std;
class Molecule;

/**
    \brief information about an atom
 
    The Atom class stores atom id, name, type, mass, charge, residue number
    
*/
class Atom
{
public:   
    /**
     * constructor
     */
    Atom(Molecule *owner, int id, string name, int resnr, double m, double q)
        : _mol(owner), _id(id), _name(name), _resnr(resnr), _m(m), _q(q)
    {
            _bPos=false;
    }
    /**
     * destructor
     */
    ~Atom() {
    }

    /**
     * get the id of the atom
     * \return bead id
     */
    const int &getId() const { return _id; }
    
    /**
     * get atom name
     * \return atom name
     */
    const string &getName() const { return _name; }
    
    /**
     * get the atom type
     * \return atom type 
     */
    const string &getType() const { return _type; }

    /**
     * get the residue number of the atom
     * \return residue id
     */
    const int &getResnr() const { return _resnr; }

    /**
     * get the charge of the atom
     * \return atom charge
     */
    const double &getQ() const { return _q; }
    
    /**
     * set the charge of the bead
     * \param q bead charge
     */
    void setQ(const double &q) { _q=q; }

    /**
     * set the position of the atom
     * \param r atom position
     */
    void setPos(const vec &r);

    /**
     * get the position of the atom
     * \return atom position
     */
    const vec &getPos() const;
      
    /**
     * direct access (read/write) to the position of the atom
     * \return reference to position 
     */
    vec &Pos() { return _pos; }

    /** does this configuration store positions? */
    bool HasPos() {return _bPos; }
           
    /** dose the bead store a position */
    void HasPos(bool b);

  
    /**
     * molecule the bead belongs to
     * \return Molecule object
     */
    Molecule *getMolecule() { return _mol; }

protected:
    int _id;
    Molecule *_mol;
    
    string _name;
    string _type;
    int _resnr;
    double _m;
    double _q;
    vec _pos;
    bool _bPos;
        
};

inline void Atom::setPos(const vec &r)
{
    _bPos=true;
    _pos = r;
}

inline const vec &Atom::getPos() const
{
    assert(_bPos);
    return _pos;
}

inline void Atom::HasPos(bool b)
{
    _bPos=b;
}

}}

#endif	/* __VOTCA_CTP_ATOM_H */

