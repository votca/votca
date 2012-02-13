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

#include <cstdlib>
#include <string>
#include <vector>

namespace votca { namespace ctp {

using namespace std;

class Topology;
class Segment;


class SegmentType
{
public:
    
    SegmentType() { }
   ~SegmentType() { _torbNrs.clear(); }

    SegmentType(int typeId,           string name,
                string basis,         string orbFile,
                string qmCoordsFile,  vector<int> torbNrs,
                bool canRigidify)
              : _id(typeId),          _name(name),
                _basisSetName(basis), _orbitalsFile(orbFile),
                _qmCoordsFile(qmCoordsFile), _torbNrs(torbNrs),
                _canRigidify(canRigidify)  { };

    SegmentType(int typeId, string name) : _id(typeId), _name(name) { };

    void setTopology(Topology *top) { _top = top; }
    void setTOrbNrs(vector<int> &torbNrs) { _torbNrs = torbNrs; }
    void setBasisName(const string &name) { _basisSetName = name; }
    void setOrbitalsFile(const string &file) { _orbitalsFile = file; }
    void setQMCoordsFile(const string &file) { _qmCoordsFile = file; }
    void setCanRigidify(bool yesno) { _canRigidify = yesno; }

    Topology     *getTopology() { return _top; }
    const string &getName() { return _name; }
    const int    &getId() { return _id; }
    vector<int>  &getTOrbNrs() { return _torbNrs; }
    const string &getBasisName() { return _basisSetName; }
    const string &getOrbitalsFile() { return _orbitalsFile; }
    const string &getQMCoordsFile() { return _qmCoordsFile; }
    const bool    canRigidify() { return _canRigidify; }



private:

    Topology*   _top;

    string      _name;
    int         _id;

    string          _qmCoordsFile;
    string          _basisSetName;
    string          _orbitalsFile;
    vector<int>     _torbNrs;

    bool            _canRigidify;


};

/*
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
    SegmentType(
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
*/



}}

#endif	/* __VOTCA_SEGMENT_TYPE_H */

