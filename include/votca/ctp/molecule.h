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

#ifndef VOTCA_CTP_MOLECULE_H
#define VOTCA_CTP_MOLECULE_H

#include <string>
#include <vector>
#include <votca/ctp/atom.h>
#include <votca/ctp/segment.h>

namespace votca { namespace ctp {
/**
    \brief Molecule is a container for atoms

    The Molecule class stores atoms. It substitutes mol_and_orbs object of libmoo.

*/
class Molecule {
public:
    /// Constructor
    Molecule(int id, string name) 
        : _id(id), _name(name)
    {}
    /// Destructor
    ~Molecule();
    /// Returns molecule ID
    const int &getId();
    /// Returns the name of the molecule
    const string &getName();



    /// Adds a pointer to a segment to this molecule
    void AddSegment( Segment* segment );    
    vector< Segment* > &Segments() { return _segments; }



    //  Add atom to internal container
    void AddAtom( Atom* atom);
    /// Returns a pointer to the atom
    Atom *getAtom(const int &id);
    /// Returns atom type
    const string &getAtomType(const int &id);
    /// Returns atom position
    const vec getAtomPosition(const int &id);
    /// Returns number of atoms in the molecule
    int NumberOfAtoms();
    /// Writes a PDB file


    void WritePDB( ostream & out );
    /// Load molecule coordinates from a file
    void ReadXYZ ( string filename );
    
    
private:
    // molecule name
    string _name ;
    // molecule ID
    int    _id;
    // map of atom names and their IDs
    map<string, Atom* > _map_AtomName_Atom;
    // list of atoms
    vector < Atom* > _atoms ;
    // list of segments
    vector < Segment* > _segments;
};

}}

#endif // VOTCA_CTP_MOLECULE_H
