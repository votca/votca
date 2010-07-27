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

#ifndef _ESPTOPOLOGYREADER_H
#define	_ESPTOPOLOGYREADER_H

#include <string>
#include "topologyreader.h"


namespace votca { namespace csg {
using namespace std;
    
/**
    \brief reader for espresso blockfiles

    This class encapsulates the espresso reading functions and provides an interface to fill a topolgy class

*/
class ESPTopologyReader
    : public TopologyReader
{
public:
    /// read a topology file
    bool ReadTopology(string file, Topology &top);
    
private:
		/// extract variable from input line
		string ReadBlockfileLine(string input_line, string variable);		
};

}}

#endif	/* _ESPTOPOLOGYREADER_H */

