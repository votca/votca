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
#include <atom.h>

namespace votca { namespace ctp {
/**
    \brief Molecule is a container for atoms

    The Molecule class stores atoms.

*/

class Molecule {
public:
    // Constructor
    Molecule() {}
    // Destructor
    ~Molecule() {}
    /// Returns molecule ID
    int getId() const {
        return _id;
    }
    /// Returns the name of the molecule
    const string &getName() const {
        return _name;
    }
    /// Returns a pointer to the atom
    Atom *getAtom(const int &id);
    /// Returns atom type
    const string &getAtomType(const int &id);
    /// Returns atom position
    const vec getAtomPosition(const int &id);
    /// Returns number of atoms in the molecule
    int NumberOfAtoms();
    /// Adds an atom the molecule
    void AddAtom(Atom *getAtom);
    /// Writes a PDB file
    void WritePDB();

private:
    // map of atom names and their IDs
    map<string, int> _names;
    // list of atoms
    vector < Atom* > _atoms ;
    // molecule name
    string _name ;
    // molecule ID
    int    _id;
};

}}

#endif // VOTCA_CTP_MOLECULE_H
