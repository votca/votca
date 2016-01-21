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

#ifndef _CRGUNITTYPE_H
#define	_CRGUNITTYPE_H

//#define DEBUG

#include "mol_and_orb.h"
#include "orbitals.h"
#include <string>
#include <stdexcept>

namespace votca { namespace moo {

class CrgUnit;

/**
    \brief information about a CrgUnitTpe

    The CrgUnitType class describes a known charge unit type. It stores
    the ID, reorganization energie, internal coordinates and orbitals.

    \todo mod orbitals to be "stripped down" for memory?
*/

class CrgUnitType
{
public:
    
    // initialised by a list of ID names, molecule specifications
    // and orbital specifications
    
   ~CrgUnitType();
       
    const vector <unsigned int>& GetTransOrbs() const;
    
    const unsigned int & getId() const { return _id; }
    
    string & GetMolName();   
    
    const string & GetName() const { return _name; }
    
    // TODO: mol_and_orb should not be in the interface!
    mol_and_orb & GetCrgUnit();

    // TODO: orb should not be in the interface!
    const orb & GetOrb() const;

    basis_set & GetBS(){ return _bs; }

    Property *getOptions() { return _options; }
    void setOptions(Property *options) { _options = options; }

private:

    basis_set               _bs;
    mol_and_orb             _intcoords;
    orb                     _orbitals;
    Property               *_options;
    vector <unsigned int>   _transorbs;
    unsigned int            _id;

    /// the name of this transport unit
    string                  _name;
    // a list of atoms in the same bead of which atoms correspond to which bead
    vector <vector <int> >  _list_atoms_monomer;

    // a list of centre of mass of the monomer
    vector < vec  >         _list_coms_monomer;

    // a list of rotation matrices to project the internal
    // coordinates of each monomer onto the reference state
    vector < matrix >       _list_ors_monomer;

    // this will take each bead and move it to positions[i] rotating by the
    // orientation corresponding to norm and rotate we assume that the pointer
    // is to a suitable molecule...
    void rotate_each_bead(vector < vec >::iterator it_pos,
                          vector < vec >::iterator it_norm,
                          vector <vec >::iterator it_plan,
                          mol_and_orb * rotated_molecule );

    // can only be created by JCalc
    CrgUnitType(): _id(-1) {};
    
    // can only be created by JCalc
    CrgUnitType(const char * namecoord,
                const char * nameorb,
                const char * nameneutr,
                const char * namecrg,
                string & basisset,
                const vector < int>& transorbs,
                const unsigned int &id,
                string name,
                vector < vector < int > > list_atoms_monomer,
                vector < vector < double > > list_weights_monomer);

    friend class CrgUnit;
    friend class JCalc;
};

inline const vector <unsigned int>& CrgUnitType::GetTransOrbs() const{
    return _transorbs;
}

inline mol_and_orb & CrgUnitType::GetCrgUnit() {
    return _intcoords;
}

inline const orb & CrgUnitType::GetOrb() const {
    return _orbitals;
}

}}

#endif	/* _CRGUNITTYPE_H */

