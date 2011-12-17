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

namespace votca { namespace ctp {

Topology::Topology()
    : _db_id(0)
{}

Topology::~Topology()
{    
    // clean up the list of fragments 
    cout << "  topology: deleting fragments: ";
    vector < Fragment* > :: iterator fragment;
    for (fragment = _fragments.begin(); fragment < _fragments.end(); ++fragment){
         cout << (*fragment)->getId() << ":" << (*fragment)->getName() << " ";      
        delete *fragment;
    }
    cout << endl;
    _fragments.clear();

    // clean up the list of segments 
    cout << "  topology: deleting segments: ";
    vector < Segment* > :: iterator segment;
    for (segment = _segments.begin(); segment < _segments.end(); ++segment){
         cout << (*segment)->getId() << ":" << (*segment)->getName() << " ";      
        delete *segment;
    }
    cout << endl;
    _segments.clear();
    
    // clean up the list of molecules; this also deletes atoms
    cout << "  topology: deleting molecules: ";
    vector < Molecule* > :: iterator molecule;
    for (molecule = _molecules.begin(); molecule < _molecules.end(); ++molecule){
         cout << (*molecule)->getId() << ":" << (*molecule)->getName() << " " ;      
        delete *molecule;
    }
    cout << endl;
    _molecules.clear();

     // clean up the list of fragment types 
    cout << "  topology: deleting fragment types: ";
    for (fragment = _fragment_types.begin(); fragment < _fragment_types.end(); ++fragment){
         cout << (*fragment)->getId() << ":" << (*fragment)->getName() << " ";      
        delete *fragment;
    }
    cout << endl;
    _fragment_types.clear();

    // clean up the list of segment types 
    cout << "  topology: deleting segment types: ";
    for (segment = _segment_types.begin(); segment < _segment_types.end(); ++segment){
         cout << (*segment)->getId() << ":" << (*segment)->getName() << " ";      
        delete *segment;
    }
    cout << endl;
    _segment_types.clear();
   
      // clean up the list of molecule types; this also deletes atoms
    cout << "  topology: deleting molecule types: ";
    for (molecule = _molecule_types.begin(); molecule < _molecule_types.end(); ++molecule){
         cout << (*molecule)->getId() << ":" << (*molecule)->getName() << " " ;      
        delete *molecule;
    }
    cout << endl;
    _molecule_types.clear();
  
    // clean up the map of the molecule md name : pointer to molecule type
    _map_MoleculeName_MoleculeType.clear();
    _map_MoleculeMDName_MoleculeName.clear();
}


Fragment *Topology::AddFragment(int fragment_id, string fragment_name, Segment* segment)
{
    Fragment* fragment = new Fragment(fragment_id, fragment_name);
    _fragments.push_back(fragment);
}

Segment *Topology::AddSegment(int segment_id, string segment_name)
{
    Segment* segment = new Segment(segment_id, segment_name);
    _segments.push_back(segment);

}

Atom *Topology::AddAtom(int atom_id, string atom_name)
{
    //_atoms.push_back( atom );
}

Molecule *Topology::AddMolecule(int molecule_id, string molecule_name)
{
    Molecule *molecule = new Molecule(molecule_id++, molecule_name);
    _molecules.push_back(molecule);
    return molecule;
}


Atom *Topology::AddAtomType(Molecule *owner, int atom_id, string atom_name, 
        int residue_number, double weight)
{
     Atom* atom = new Atom(owner, atom_id, atom_name, residue_number, weight);
     _atom_types.push_back(atom);  
    return atom;
}

Fragment *Topology::AddFragmentType(int fragment_id, Property *property)
{
    string fragment_name = property->get("name").as<string>();    
    Fragment* fragment = new Fragment(fragment_id, fragment_name);
     _fragment_types.push_back(fragment);  
    return fragment;
}

Segment *Topology::AddSegmentType(int segment_id, Property *property)
{
    string segment_name = property->get("name").as<string>();
    Segment* segment = new Segment(segment_id, segment_name);
    _segment_types.push_back(segment);
    return segment;

}

Molecule *Topology::AddMoleculeType(int molecule_id, Property *property)
{
    string molecule_name = property->get("name").as<string>();
    string molecule_mdname = property->get("mdname").as<string>();

    Molecule *molecule = new Molecule(molecule_id, molecule_name);
    _molecule_types.push_back(molecule);
    // TODO check if the name is already there
    _map_MoleculeMDName_MoleculeName[molecule_name] = molecule_mdname;
    _map_MoleculeName_MoleculeType[molecule_name] = molecule;
    return molecule;
}

Molecule *Topology::getMoleculeType(string name)
{
    return _map_MoleculeName_MoleculeType.at(name);
}
   
   
void Topology::ParseSegmentDefinitions( Property &topology )
{

    string fragment_name;
   
    if ( tools::globals::verbose ) {
        cout << "Topology: Parsing the partitioning on segments and fragments" << endl;
        //cout << " molecule[name:id] segment[name:type] fragment[name:type:position] " << endl;
    }
    
    
    list<Property *> molecules = topology.Select("topology.molecules.molecule"); 
    list<Property *>::iterator it_molecule;

    cout << " Found " << molecules.size() << " molecule types" << endl;
            
    // load all coordinates of atoms in molecules, create fragments and segments
    int molecule_id = 1;
    for ( it_molecule = molecules.begin(); it_molecule != molecules.end(); ++it_molecule ){
        
        Molecule *molecule = AddMoleculeType(molecule_id++, *it_molecule );
        
        // load the segments
        list<Property *> segments = (*it_molecule)->Select("segments.segment"); 
        list<Property *>::iterator it_segment;
        cout << "  - Found " << segments.size() << " segments in this molecule type" << endl;
        
        int segment_id = 1;
        for ( it_segment = segments.begin(); it_segment != segments.end(); ++it_segment ){
                // Create a new segment
                Segment *segment = AddSegmentType(segment_id++, *it_segment );
                molecule->AddSegment( segment );
                
                // load the fragments
                list<Property *> fragments = (*it_segment)->Select("fragments.fragment"); 
                list<Property *>::iterator it_fragment;
                cout << "    - Found " << fragments.size() << " fragments in this segment" << endl;
                
                int fragment_id = 1;
                for ( it_fragment = fragments.begin(); it_fragment != fragments.end(); ++it_fragment ){
                    Fragment* fragment = AddFragmentType(fragment_id, *it_fragment);
                    segment->AddFragment( fragment );
                    
                    string mdatoms = (*it_fragment)->get("mdatoms").as<string>();
                    string qmatoms = (*it_fragment)->get("qmatoms").as<string>();
                    string weights = (*it_fragment)->get("weights").as<string>();
                    
                    Atom *atom = AddAtomType(molecule, atom_id, atom_name, residue_number, weight);
                    fragment->AddAtom( atom );
                    segment->AddAtom( atom );
                    

                    
                }
        }
       
    }
    
    /*
    if ( tools::globals::verbose ) {
        cout << " segment[" << segment_name  << ":" << segment_type  << "] "
             << "fragment[" << bead_name  << ":" << bead_type << ":" << bead_position  << "] "
            << "molecule[" << molecule_name  << ":" << molecule_id  << "] " << endl;
    }
     */
 
    if ( tools::globals::verbose ) {
        //cout << " segment[name:type] fragment[name:type:position] molecule[name:id] " << endl;
        cout << "Topology: Done with parsing the partitioning" << endl;
    }

}


}}
