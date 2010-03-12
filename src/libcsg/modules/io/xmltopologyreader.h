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

#ifndef _XMLTOPOLOGYREADER_H
#define	_XMLTOPOLOGYREADER_H

#include <string>
#include "topologyreader.h"
#include <stack>

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
    void ParseIgnore(const string &el, map<string, string> &attr);
    
public:


    /// start element callback for xml parser
    void StartElemHndl(const string &el, map<string, string> &attr);
    /// end element callback for xml parser
    void EndElemHndl(const string &el);

private:

    Topology *_top;

    typedef void (XMLTopologyReader::*ElemHandler_t)(const string &, map<string, string> &);

    void SetHandler(ElemHandler_t handler);

    stack<ElemHandler_t> _stack_handler;
    ElemHandler_t _handler;
};

#endif	/* _PDBTOPOLOGYREADER_H */

