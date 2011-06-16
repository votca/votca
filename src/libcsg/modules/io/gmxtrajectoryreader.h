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

#ifndef _gmxtrajectoryreader_H
#define	_gmxtrajectoryreader_H

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include <string>
#include <votca/csg/trajectoryreader.h>
#include "gmx_version_check.h"

#if GMX == 50
        #include <gromacs/legacyheaders/statutil.h>
        #include <gromacs/legacyheaders/typedefs.h>
        #include <gromacs/legacyheaders/smalloc.h>
        #include <gromacs/legacyheaders/vec.h>
        #include <gromacs/legacyheaders/copyrite.h>
        #include <gromacs/legacyheaders/statutil.h>
        #include <gromacs/legacyheaders/tpxio.h>
#elif GMX == 45
        #include <gromacs/statutil.h>
        #include <gromacs/typedefs.h>
        #include <gromacs/smalloc.h>
        #include <gromacs/vec.h>
        #include <gromacs/copyrite.h>
        #include <gromacs/statutil.h>
        #include <gromacs/tpxio.h>
#elif GMX == 40
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
#else
#error Unsupported GMX version
#endif
    // this one is needed because of bool is defined in one of the headers included by gmx
    #undef bool

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;

/**
    \brief class for reading gromacs trajectory files

    This class provides the TrajectoryReader interface and encapsulates the trajectory reading function of gromacs

*/
class GMXTrajectoryReader : public TrajectoryReader
{
    public:
        GMXTrajectoryReader() {
            gmx::CheckVersion();
        }

        /// open a trejectory file
        bool Open(const string &file);
        /// read in the first frame
        bool FirstFrame(Topology &top);
        /// read in the next frame
        bool NextFrame(Topology &top);

        void Close();
        
    private:
        string _filename;
        
        // gmx status used in read_first_frame and _read_next_frame;
#if GMX == 50
       t_trxstatus* _gmx_status;
#elif GMX == 45
       t_trxstatus* _gmx_status;
#elif GMX == 40
       int _gmx_status;
#else
#error Unsupported GMX version
#endif
        /// gmx frame
        t_trxframe _gmx_frame;
        
};

}}

#endif	/* _gmxtrajectoryreader_H */

