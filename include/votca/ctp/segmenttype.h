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

#ifndef __VOTCA_SEGMENT_TYPE_H
#define	__VOTCA_SEGMENT_TYPE_H

#include <string>
#include <votca/ctp/segment.h>

namespace votca { namespace ctp {

/**
    \brief Type of a conjugated segment

    Stores the ID, reorganization energy, and internal coordinates.

*/
class SegmentType
{
public:
    /// Constructor
    SegmentType() {};
    /// Destructor  
    ~SegmentType() {};
    /// Segment id
    const unsigned int & getId() const { return _id; }
    /// Name of the molecule this segment belongs to
    const string & getMoleculeName() {};   
    /// Name of the conjugated segment
    const string & getName() const { return _name; }
    /// Options provided in the xml file
    Property    *getOptions() { return _options; }
    /// Change options provided in the xml file
    void        setOptions(Property *options) { _options = options; }

private:

    Property *_options;
    /// ID of the conjugated segment
    unsigned int      _id;
    /// Name of the conjugated segment
    string            _name;  
    // List of atoms 
    vector < vector <int>  > _list_atoms_monomer;
    /// Coordinates of fragments 
    vector < vec  >   _list_coms_monomer;
    /// List of rotation matrices to put the internal coordinates of each 
    /// monomer onto the reference state
    vector < matrix > _list_ors_monomer;

    /// this willA take each bead and move it to positions[i] rotating by the
    /// orientation corresponding to norm and rotate we assume that the pointer
    /// is to a suitable molecule...
    void SegmentType(
            vector < vec >::iterator it_pos , vector < vec >::iterator it_norm,
            vector <vec >::iterator it_plan, Segment *segment );

    
    /// Constructor
    SegmentType(
            const char * namecoord, const char * nameorb, const char * nameneutr,
            const char * namecrg, string & basisset, 
            const vector < int>& transorbs,
            const unsigned int &id,  string name,
            vector < vector < int > > list_atoms_monomer,
            vector < vector < double > > list_weights_monomer);

    friend class Segment;
};

}}

#endif	/* __VOTCA_SEGMENT_TYPE_H */

