/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_CSG_TRAJECTORYREADER_H
#define	_VOTCA_CSG_TRAJECTORYREADER_H

#include <string>
#include "topology.h"
#include "fileformatfactory.h"

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;

/**
    \brief trajectoryreader interface
    
    This class defines the interface a trajectory reader has to implement
 */
class TrajectoryReader
{
public:
    virtual ~TrajectoryReader() {}
    /// open a trejectory file
    virtual bool Open(const string &file) = 0;
        
    virtual void Close() {};
        
    /// read in the first frame
    virtual bool FirstFrame(Topology &top) = 0;
    /// read in the next frame
    virtual bool NextFrame(Topology &top) = 0;

    static void RegisterPlugins(void);
};

// important - singleton pattern, make sure factory is created before accessed
inline FileFormatFactory<TrajectoryReader> &TrjReaderFactory()
{
    static FileFormatFactory<TrajectoryReader> _TrjReaderFactory;
    return _TrjReaderFactory;
}

}}

#endif	/* _VOTCA_CSG_TRAJECTORYREADER_H */

