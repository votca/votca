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

#ifndef __VOTCA_CSG_XYZREADER_H
#define	__VOTCA_CSG_XYZREADER_H

#include <string>
#include <iostream>
#include <fstream>
#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>

namespace votca { namespace csg {

/**
    \brief class for reading xyz files

    This class provides the TrajectoryReader + Topology reader interface
    for xyz files

*/
class XYZReader :
    public TrajectoryReader, public TopologyReader
{
    public:
        XYZReader() {}
        ~XYZReader() {}

        /// open a topology file
        bool ReadTopology(std::string file, Topology &top);

        /// open a trejectory file
        bool Open(const std::string &file);
        /// read in the first frame
        bool FirstFrame(Topology &top);
        /// read in the next frame
        bool NextFrame(Topology &top);

        void Close();

    private:

        template <bool topology>
        bool ReadFrame(Topology &top);

        std::ifstream _fl;

        int _line;
};

}}

#endif

