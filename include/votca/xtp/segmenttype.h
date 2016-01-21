/*
 *            Copyright 2009-2016 The VOTCA Development Team
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


#ifndef __VOTCA_SEGMENT_TYPE_H
#define	__VOTCA_SEGMENT_TYPE_H

#include <cstdlib>
#include <string>
#include <vector>

namespace votca { namespace xtp {

using namespace std;

class Topology;
class Segment;


class SegmentType
{
public:
    
    SegmentType() { }
   ~SegmentType() { }

    SegmentType(int typeId,                  string name,
                string basis,                string orbFile,
                string qmCoordsFile,         bool canRigidify)
              : _id(typeId),                 _name(name),
                _basisSetName(basis),        _orbitalsFile(orbFile),
                _qmCoordsFile(qmCoordsFile), _canRigidify(canRigidify)
              { };

    SegmentType(int typeId, string name) : _id(typeId), _name(name) { };

    void setTopology(Topology *top) { _top = top; }
    void setBasisName(const string &name) { _basisSetName = name; }
    void setOrbitalsFile(const string &file) { _orbitalsFile = file; }
    void setQMCoordsFile(const string &file) { _qmCoordsFile = file; }
    void setCanRigidify(bool yesno) { _canRigidify = yesno; }

    Topology     *getTopology() { return _top; }
    const string &getName() { return _name; }
    const int    &getId() { return _id; }
    const string &getBasisName() { return _basisSetName; }
    const string &getOrbitalsFile() { return _orbitalsFile; }
    const string &getQMCoordsFile() { return _qmCoordsFile; }
    const bool   &canRigidify() { return _canRigidify; }



private:

    Topology*   _top;

    int         _id;
    string      _name;

    string          _basisSetName;
    string          _orbitalsFile;
    string          _qmCoordsFile;

    bool            _canRigidify;


};

}}

#endif	/* __VOTCA_SEGMENT_TYPE_H */

