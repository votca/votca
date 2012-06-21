/*
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef VOTCA_CTP_MOLECULE_H
#define VOTCA_CTP_MOLECULE_H

#include <string>
#include <vector>
#include <votca/ctp/atom.h>
#include <votca/ctp/segment.h>
#include <votca/ctp/fragment.h>

class Topology;

namespace votca { namespace ctp {

class Molecule {
public:

    Molecule(int id, string name) : _id(id), _name(name) {}
    Molecule() { }
   ~Molecule();

    const int       &getId();
    const string    &getName();

    void AddSegment( Segment* segment );
    void AddFragment( Fragment* fragment);
    void AddAtom( Atom* atom);

    vector< Atom* >     &Atoms() { return _atoms; }
    vector< Fragment* > &Fragments() { return _fragments; }
    vector< Segment* >  &Segments() { return _segments; }

    Atom           *getAtom(const int &id);
    const string   &getAtomType(const int &id);
    const vec       getAtomPosition(const int &id);
    int             NumberOfAtoms();

    inline void setTopology(Topology *container) { _top = container; }
    Topology   *getTopology() { return _top; }

    /// Load molecule coordinates from a file
    void ReadXYZ ( string filename );
    void WritePDB( FILE *out );
    
private:

    Topology *_top;

    vector < Segment* >   _segments;
    vector < Fragment* >  _fragments;
    vector < Atom* >      _atoms ;

    string  _name ;
    int     _id;

    map<string, Atom* > _map_AtomName_Atom;

};

}}

#endif // VOTCA_CTP_MOLECULE_H
