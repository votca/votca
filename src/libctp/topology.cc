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
    
}


void Topology::Initialize( Property &topology )
{

    string coordinate_file;
    string molecule_name;
    string segment_name;
    string fragment_name;
   
    if ( tools::globals::verbose ) {
        cout << "Initializing QM molecules" << endl;
        cout << " molecule[name:id] segment[name:type] fragment[name:type:position] " << endl;
    }
    
    
    list<Property *> molecules = topology.Select("topology.molecules.molecule"); 
    list<Property *>::iterator it_molecule;

    cout << " Found " << molecules.size() << " molecule types" << endl;
            
    // load all coordinates of atoms in molecules, create fragments and segments
    int molecule_id = 1;
    for ( it_molecule = molecules.begin(); it_molecule != molecules.end(); ++it_molecule ){
        
        molecule_name = (*it_molecule)->get("name").as<string>();
        Molecule *molecule = new Molecule(molecule_id++, molecule_name);
        _molecules.push_back(molecule);        
        
        // load the segments
        list<Property *> segments = (*it_molecule)->Select("segments.segment"); 
        list<Property *>::iterator it_segment;
        cout << "  - Found " << segments.size() << " segments in this molecule" << endl;
        
        int segment_id = 1;
        for ( it_segment = segments.begin(); it_segment != segments.end(); ++it_segment ){
                // Create a new segment
                segment_name = (*it_segment)->get("name").as<string>();
                Segment* segment = new Segment( segment_id++, segment_name );
                _segments.push_back(segment);

                // load the fragments
                list<Property *> fragments = (*it_segment)->Select("fragments.fragment"); 
                list<Property *>::iterator it_fragment;
                cout << "    - Found " << fragments.size() << " fragments in this segment" << endl;
                
                int fragment_id = 1;
                for ( it_fragment = fragments.begin(); it_fragment != fragments.end(); ++it_fragment ){
                    fragment_name = (*it_fragment)->get("name").as<string>();
                    Fragment* fragment = new Fragment( fragment_id++, fragment_name );
                    _fragments.push_back(fragment);

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
        cout << " segment[name:type] fragment[name:type:position] molecule[name:id] " << endl;
        cout << "Done with initializing conjugated segments" << endl;
    }

}


}}
