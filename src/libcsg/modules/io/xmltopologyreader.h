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
// File:   pdbtopologyreader.h
// Author: victor
//
// Created on 14. April 2008, 13:33
//

#ifndef _XMLTOPOLOGYREADER_H
#define	_XMLTOPOLOGYREADER_H

#include <string>
#include <libxml/xmlreader.h>
#include "topologyreader.h"
#include <stack>

using namespace std;
    
/**

*/
class XMLTopologyReader
   : public TopologyReader
{
public:
    /// read a topology file
    bool ReadTopology(string file, Topology &top);

private:    
    void ReadTopolFile(string file);
          
    void ParseTopology(xmlNodePtr node);
    void ParseMolecules(xmlNodePtr node);
    
    Topology *_top;
};

#endif	/* _PDBTOPOLOGYREADER_H */

