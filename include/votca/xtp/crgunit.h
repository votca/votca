/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _CRGUNIT_H
#define	_CRGUNIT_H

#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>
#include "global.h"
#include "crgunittype.h"
#include "units.h"

namespace votca { namespace xtp {

using namespace std;

/**
 * \brief information about a CrgUnit
 *
 * The CrgUnit class describes a charge unit. It stores information like
 * the id, the position, orientation, type of the charge transport unit.
 *
 */

class CrgUnit
{
public:

    CrgUnit() : _id(-1), _type(NULL), _molid(-1)  {};

    virtual ~CrgUnit();

    CrgUnit(vector <vec> positions,
            vector <vec> norms,
            vector <vec> planes,
            const unsigned int & id, CrgUnitType * type,
            const unsigned int & molId);

    CrgUnit(const unsigned int & id,
            CrgUnitType * type,
            const unsigned int & molId);
    
    /// TODO: vr have a lock
    void copyCrgUnit(CrgUnit & acrg);
    /// TODO: vr have a lock
    void copyCrgUnit(CrgUnit & acrg, const int & id);    

    const unsigned int& getId() const;
    void setId(const int& i);

    CrgUnitType* getType() const;

    const string &getName() const { return _name; }
    void setName(const string &name) { _name = name; }

    const unsigned int& getMolId() const;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    int GetN();

    void SetPos(const unsigned int& i, const vec& pos);
    void SetNorm(const unsigned int& i, const vec& pos);
    void SetPlane(const unsigned int& i, const  vec& pos);

    vec GetPos(const int & i);
    vec GetNorm(const int & i);
    vec GetPlane(const int & i);
    vec GetCom();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // what of this can we simplify/beautify?
    mol_and_orb * rotate_translate_beads();

    // this is the function called from the rate
    // calculator on already initialised molecules
    void rot_two_mol(CrgUnit & two, mol_and_orb & mol1, mol_and_orb & mol2);
    void shift(const vec & displ);
    void rotate(matrix mat);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

private:
    // variables are internally stored in atomic units (Hartree, Bohr etc.),
    // but printouts are always nm, eV, SI.
    unsigned int    _id;
    // the type
    string          _name;
    // vector of CoMs of monomers (stored in Bohr, returned in nm)
    vector < vec >  _positions;
    // normal vector of monomers (normalized)
    vector < vec >  _norms;
    // orientation vector of monomers (normalized)
    vector < vec >  _planes;
    // a reference to the crgunittype
    CrgUnitType *   _type;
    // the molecule index
    unsigned int    _molid;

    vector <vec> shift_pos(const vec & a);
};

inline CrgUnit::CrgUnit(const unsigned int & id,
                        CrgUnitType * type,
                        const unsigned int & molId)
                      : _id(id), _type(type), _molid(molId) { }

inline void CrgUnit::copyCrgUnit(CrgUnit & acrg, const int & id) {
    copyCrgUnit(acrg);
    _id = id;
}

inline const unsigned int& CrgUnit::getId() const { return _id; }
inline void CrgUnit::setId(const int& i) { _id = i; }

inline CrgUnitType* CrgUnit::getType() const { return _type; }

inline const unsigned int& CrgUnit::getMolId() const { return _molid; }

inline int CrgUnit::GetN() {
    if (_norms.size() != _planes.size() ||
        _norms.size() != _positions.size()) {
        cerr << "Error in the crg unit of id " << _id << endl;
    }
    return _norms.size();
}

inline void CrgUnit::SetPos(const unsigned int& i, const vec& pos) {
    if ( i >= _positions.size() ) {
         _positions.resize(i+1);
         _norms.resize(i+1);
         _planes.resize(i+1);
    }
    _positions[i] = unit<nm,bohr>::to(pos);
};

inline void CrgUnit::SetNorm(const unsigned int& i, const vec& pos) {
    if ( i >= _norms.size() ) {
         _positions.resize(i+1);
         _norms.resize(i+1);
         _planes.resize(i+1);
    }
    _norms[i] = pos;
};

inline void CrgUnit::SetPlane(const unsigned int& i, const vec& pos) {
    if ( i >= _planes.size() ) {
         _positions.resize(i+1);
         _norms.resize(i+1);
         _planes.resize(i+1);
    }
    _planes[i] = pos;
};

inline vec CrgUnit::GetPos(const int & i) {
    return unit<bohr,nm>::to(_positions[i]);
}

inline vec CrgUnit::GetNorm(const int & i) {
    return _norms[i];
}

inline vec CrgUnit::GetPlane(const int & i) {
    return _planes[i];
}

inline  vec  CrgUnit::GetCom()  {
    vector <vec>::iterator itp = _positions.begin();
    vec com(0.,0.,0.);
    for (  ; itp!=_positions.end() ; ++itp )
           com += *itp;
    return unit<bohr,nm>::to( com/double(_positions.size()) );
}

inline void CrgUnit::shift(const vec & displ) {
    vec displ_bohr = unit<nm,bohr>::to(displ);
    vector <vec> ::iterator it_vec;
    for (it_vec = _positions.begin();
         it_vec != _positions.end();
         ++it_vec)
         *it_vec = *it_vec + displ_bohr;
}

inline vector <vec> CrgUnit::shift_pos(const vec & a) {
    vector <vec>::iterator it_pos;
    vector <vec> res;
    for (it_pos = _positions.begin(); 
         it_pos != _positions.end();
         ++it_pos) {
         vec b = (*it_pos) + a;
         res.push_back(b);
    }
    return res;
}

}}

#endif	/* _CRGUNIT_H */

