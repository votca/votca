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

#ifndef __VOTCA_CTP_TOPOLOGY_H
#define	__VOTCA_CTP_TOPOLOGY_H

#include <votca/tools/property.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/fragment.h>
#include <votca/ctp/segment.h>
#include <votca/ctp/molecule.h>

namespace votca { namespace ctp {

/**
 * \brief Container for molecules, conjugated segments, rigid fragments,
 * and atoms.
*/
class Topology 
{
public:
    Topology();
   ~Topology();

    /// Load the topology based on definitions of conjugated segments
    void ParseSegmentDefinitions( Property &topology );

     /// Creates a fragment and adds it to the topology
    Fragment *AddFragment (int fragment_id, string fragment_name, Segment* segment);
    /// Creates a segment and adds it to the topology
    Segment *AddSegment (int segment_id, string segment_name);   
    /// Creates an atom and adds it to the topology
    Atom *AddAtom (int atom_id, string atom_name);   
    /// Creates a molecule and adds it to the topology
    Molecule *AddMolecule (int molecule_id, string molecule_name);
   
    int getDatabaseId() { return _db_id; };
    void setDatabaseId(int id) { _db_id = id; }

   
protected:

    vector < Molecule* > _molecules;
    vector < Segment* > _segments;
    vector < Fragment* > _fragments;
    vector < Atom* > _atoms;
    
    vector < Molecule* > _molecule_types;
    vector < Segment* > _segment_types;
    vector < Fragment* > _fragment_types;
    vector < Atom* > _atom_types;
    
    map < string, Molecule* > _map_MoleculeName_MoleculeType;
    map < string, string > _map_MoleculeMDName_MoleculeName;
    
    map < int, Segment* > _map_id_segment;


    int _db_id;
    /// Adds an atom type when parsing segment/fragment definitions
    Atom *AddAtomType(Molecule *owner, int atom_id, string atom_name, 
        int residue_number, double weight);   
    /// Adds a fragment type for internal use (when loading segments.xml) 
    Fragment *AddFragmentType(int fragment_id, Property *property);
    /// Adds a segment type (when loading segments.xml)
    Segment *AddSegmentType(int segment_id, Property *property);
    /// Adds a molecule type (when loading segments.xml)
    Molecule *AddMoleculeType (int molecule_id, Property *property);
    /// Returns a pointer to a molecule type with a specified name 
    Molecule *getMoleculeType(string name);
   
};

}}

#endif	/* __VOTCA_CTP_TOPOLOGY_H */

