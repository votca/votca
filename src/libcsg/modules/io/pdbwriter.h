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

#ifndef _PDBWRITER_H
#define	_PDBWRITER_H

#include <stdio.h>
#include <votca/csg/topology.h>
#include <votca/csg/trajectorywriter.h>

namespace votca { namespace csg {

class PDBWriter
: public TrajectoryWriter
{
public:
    
    void Open(std::string file, bool bAppend = false);
    void Close();
    
    void RegisteredAt(ObjectFactory<std::string, TrajectoryWriter> &factory) {}    

    void Write(Topology *conf);

private:
    FILE *_out;
};

}}

#endif	/* _PDBWRITER_H */

