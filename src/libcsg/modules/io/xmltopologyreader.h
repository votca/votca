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

#ifndef _XMLTOPOLOGYREADER_H
#define	_XMLTOPOLOGYREADER_H

#include <string>
#include "topologyreader.h"
#include <stack>
#include <votca/tools/parsexml.h>


namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;

/**
 *  Reads in an xml topology
 *
 * \todo this is a sloppy implementation using expat, is just reads attributes
 * \todo should be extended to also read beads, ...
 * 
*/
class XMLTopologyReader
   : public TopologyReader
{
public:
    /// read a topology file
    bool ReadTopology(string file, Topology &top);

private:    

    void ReadTopolFile(string file);

    void ParseRoot(const string &el, map<string, string> &attr);
    void ParseTopology(const string &el, map<string, string> &attr);
    void ParseMolecules(const string &el, map<string, string> &attr);
    
private:
    ParseXML _parser;

    Topology *_top;
};

}}

#endif	/* _PDBTOPOLOGYREADER_H */

