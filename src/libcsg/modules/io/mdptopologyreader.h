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
 * File:   mdptopologyreader.h
 * Author: ruehle
 *
 * Created on June 12, 2008, 10:47 AM
 */

#ifndef _MDPTOPOLOGYREADER_H
#define	_MDPTOPOLOGYREADER_H

#include "topologyreader.h"

/**
    \brief reader for Alexander Lyubartsev's md format
 
 */
class MDPTopologyReader
    : public TopologyReader
{
public:
    /// read a topology file
    bool ReadTopology(string file, Topology &top);
    TopologyReader *Clone() { return new MDPTopologyReader(); }
    
private:
       
};


#endif	/* _MDPTOPOLOGYREADER_H */

