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

#include <votca/ctp/topology.h>
#include <votca/tools/globals.h>
#include <boost/lexical_cast.hpp>

namespace votca { namespace ctp {

Topology::Topology() : _db_id(0) { }

Topology::~Topology() {
    // clean up the list of fragments 
    vector < Fragment* > :: iterator fragment;
    for (fragment = _fragments.begin();
         fragment < _fragments.end();
         ++fragment) {
         delete *fragment;
    }
    _fragments.clear();

    // clean up the list of segments 
    vector < Segment* > :: iterator segment;
    for (segment = _segments.begin();
         segment < _segments.end();
         ++segment) {
         delete *segment;
    }
    _segments.clear();
    
    // clean up the list of molecules; this also deletes atoms
    vector < Molecule* > :: iterator molecule;
    for (molecule = _molecules.begin();
         molecule < _molecules.end();
         ++molecule) {
         delete *molecule;
    }
    _molecules.clear();   
}


Fragment *Topology::AddFragment(int fragment_id, string fragment_name,
                                Segment* segment) {
    Fragment* fragment = new Fragment(fragment_id, fragment_name);
    _fragments.push_back(fragment);
    return fragment;
}

Segment *Topology::AddSegment(int segment_id, string segment_name) {
    Segment* segment = new Segment(segment_id, segment_name);
    _segments.push_back(segment);
    return segment;
}

Atom *Topology::AddAtom(int atom_id, string atom_name) {
    Atom *atom = new Atom;
    _atoms.push_back(atom);
    return atom;
}

Molecule *Topology::AddMolecule(int molecule_id, string molecule_name) {
    Molecule *molecule = new Molecule(molecule_id++, molecule_name);
    _molecules.push_back(molecule);
    return molecule;
}

}}
