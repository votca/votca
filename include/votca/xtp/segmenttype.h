/*
 *            Copyright 2009-2018 The VOTCA Development Team
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
/// For earlier commit history see ctp commit 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef VOTCA_XTP_SEGMENT_TYPE_H
#define	VOTCA_XTP_SEGMENT_TYPE_H

#include <string>

namespace votca { namespace xtp {

class Topology;
class Segment;

class SegmentType
{
public:
    
    SegmentType() { }
   ~SegmentType() { }

    SegmentType(int typeId,                  std::string name,
                std::string basis,                std::string orbFile,
                std::string qmCoordsFile,         bool canRigidify)
              : _id(typeId),                 _name(name),
                _basisSetName(basis),        _orbitalsFile(orbFile),
                _qmCoordsFile(qmCoordsFile), _canRigidify(canRigidify)
              { };

    SegmentType(int typeId, std::string name) : _id(typeId), _name(name) { };

    void setTopology(Topology *top) { _top = top; }
    void setBasisName(const std::string &name) { _basisSetName = name; }
    void setOrbitalsFile(const std::string &file) { _orbitalsFile = file; }
    void setQMCoordsFile(const std::string &file) { _qmCoordsFile = file; }
    void setCanRigidify(bool yesno) { _canRigidify = yesno; }

    Topology     *getTopology() { return _top; }
    const std::string &getName() { return _name; }
    const int    &getId() { return _id; }
    const std::string &getBasisName() { return _basisSetName; }
    const std::string &getOrbitalsFile() { return _orbitalsFile; }
    const std::string &getQMCoordsFile() { return _qmCoordsFile; }
    const bool   &canRigidify() { return _canRigidify; }



private:

    Topology*   _top;

    int         _id;
    std::string      _name;

    std::string          _basisSetName;
    std::string          _orbitalsFile;
    std::string          _qmCoordsFile;

    bool            _canRigidify;


};

}}

#endif	// VOTCA_XTP_SEGMENT_TYPE_H 

