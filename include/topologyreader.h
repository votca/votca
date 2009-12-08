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
// 
// File:   topologyreader.h
// Author: ruehle
//
// Created on January 18, 2008, 3:36 PM
//

#ifndef _TOPOLOGYREADER_H
#define	_TOPOLOGYREADER_H

#include <string>
#include "topology.h"
#include "fileformatfactory.h"

using namespace std;

class TopologyReader
{
public:
    virtual ~TopologyReader() {}
    /// open a trejectory file
    virtual bool ReadTopology(string file, Topology &top) = 0;
        
    static void RegisterPlugins(void);
};

// important - singleton pattern, make sure factory is created before accessed
inline FileFormatFactory<TopologyReader> &TopReaderFactory()
{
    static FileFormatFactory<TopologyReader> _TopReaderFactory;
    return _TopReaderFactory;
}

#endif	/* _TOPOLOGYREADER_H */

