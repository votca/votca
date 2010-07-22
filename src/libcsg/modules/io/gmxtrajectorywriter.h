/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _GMXTRAJECTORYWRITER_H
#define	_GMXTRAJECTORYWRITER_H

#include "topology.h"
#include "trajectorywriter.h"

#include "config.h"

#ifdef GMX4DEV
        #include <gromacs/statutil.h>
        #include <gromacs/typedefs.h>
        #include <gromacs/smalloc.h>
        #include <gromacs/vec.h>
        #include <gromacs/copyrite.h>
        #include <gromacs/statutil.h>
        #include <gromacs/tpxio.h>
#else
   extern "C"
   {
        #include <statutil.h>
        #include <typedefs.h>
        #include <smalloc.h>
        #include <vec.h>
        #include <copyrite.h>
        #include <statutil.h>
        #include <tpxio.h>
    }
#endif
    // this one is needed because of bool is defined in one of the headers included by gmx
    #undef bool

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;

class GMXTrajectoryWriter
    : public TrajectoryWriter
{
public:
    
    void Open(string file, bool bAppend = false);
    void Close();
    void Write(Topology *conf);

    private:
#ifdef GMX4DEV
       t_trxstatus* _file;
#else
       int _file;
#endif
};

}}

#endif	/* _GMXTRAJECTORYWRITER_H */

