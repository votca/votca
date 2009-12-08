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
/* 
 * File:   mdptrajectoryreader.h
 * Author: victorr
 *
 * Created on June 12, 2008, 1:51 PM
 */

#ifndef _MDPTRAJECTORYREADER_H
#define	_MDPTRAJECTORYREADER_H

#include <stdio.h>
#include "trajectoryreader.h"

class MDPTrajectoryReader : public TrajectoryReader
{
    public:        
        /// open a trejectory file
        bool Open(const string &file);
        /// read in the first frame
        bool FirstFrame(Configuration &conf);
        /// read in the next frame
        bool NextFrame(Configuration &conf);

        void Close();
        
        TrajectoryReader *Clone() { return dynamic_cast<TrajectoryReader*>(new MDPTrajectoryReader()); }

    private:
        FILE *_fl;
        int _moltypes;
        vector<int> _nmols;
        vector<int> _natoms;
};



#endif	/* _MDPTRAJECTORYREADER_H */

