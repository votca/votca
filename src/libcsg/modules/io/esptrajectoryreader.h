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

#ifndef _ESPTRAJECTORYREADER_H
#define	_ESPTRAJECTORYREADER_H

#include "trajectoryreader.h"

namespace votca { namespace csg {

using namespace std;

class ESPTrajectoryReader : public TrajectoryReader
{
    public:        
        /// open a trajectory file
        bool Open(const string &file);
        /// read in the first frame
        bool FirstFrame(Topology &top);
        /// read in the next frame
        bool NextFrame(Topology &top);

        void Close();
        
        TrajectoryReader *Clone() { return dynamic_cast<TrajectoryReader*>(new ESPTrajectoryReader()); }

    private:
        string _fl;
				string _temp_nextframe;				
};

}}

#endif	/* _ESPTRAJECTORYREADER_H */

