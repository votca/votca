/* 
 * Copyright 2009-2015 The VOTCA Development Team (http://www.votca.org)
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

#ifndef HAVE_NO_CONFIG
#include <votca_config.h>
#endif

#include <votca/csg/topology.h>
#include <votca/csg/trajectorywriter.h>

#include <gromacs/trajectory/trajectoryframe.h>
#include <gromacs/fileio/trxio.h>
// this one is needed because of bool is defined in one of the headers included by gmx
#undef bool

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;

class GMXTrajectoryWriter
    : public TrajectoryWriter
{
public:
    GMXTrajectoryWriter() {}

    void Open(string file, bool bAppend = false);
    void Close();
    void Write(Topology *conf);

    private:
       t_trxstatus* _file;
};

}}

#endif	/* _GMXTRAJECTORYWRITER_H */

