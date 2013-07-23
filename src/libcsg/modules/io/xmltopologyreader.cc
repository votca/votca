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

#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <stdio.h>
#include <stdexcept>
#include "xmltopologyreader.h"

namespace votca { namespace csg {

bool XMLTopologyReader::ReadTopology(string filename, Topology &top)
{ 
  _top = &top;

  _parser.NextHandler(this, &XMLTopologyReader::ParseRoot);
  _parser.Open(filename);
    
  return true;
}

void XMLTopologyReader::ReadTopolFile(string file)
{
    TopologyReader *reader;
    reader = TopReaderFactory().Create(file);
    if(!reader)
        throw runtime_error(file + ": unknown topology format");
    
    reader->ReadTopology(file, *_top);
    
    delete reader;
}

void XMLTopologyReader::ParseRoot(const string &el, map<string, string> &attr)
{
    if(el == "topology") {
        if(attr["base"] != "")
            ReadTopolFile(attr["base"]);
        _parser.NextHandler(this, &XMLTopologyReader::ParseTopology);
    }
    else {
        throw std::runtime_error("wrong root node in xml topology file");
    }
}

void XMLTopologyReader::ParseTopology(const string &el, map<string, string> &attr)
{
    if(el == "molecules") {
         _parser.NextHandler(this, &XMLTopologyReader::ParseMolecules);
    }
    else if(el == "box") {
        matrix m;
        m.ZeroMatrix();
        m[0][0] = boost::lexical_cast<double>(attr["xx"]);
        m[1][1] = boost::lexical_cast<double>(attr["yy"]);
        m[2][2] = boost::lexical_cast<double>(attr["zz"]);
        _top->setBox(m);
        _parser.NextHandler(this, &XMLTopologyReader::ParseTopology);
    }
    else if(el == "beadtypes") {
         _parser.NextHandler(this, &XMLTopologyReader::ParseBeadTypes);
    }
    else {
        throw runtime_error("unknown tag: "+ el);
    }
}

void XMLTopologyReader::ParseMolecules(const string &el, map<string, string> &attr)
{
    map<string, string>::iterator iter;
    if (el == "clear") {
        _top->ClearMoleculeList();
        _parser.IgnoreChilds();
    }
    else if (el == "rename") {
        string molname = attr["name"];
        string range = attr["range"];
        if (molname == "" || range == "")
            throw runtime_error("invalid rename tag");
        _top->RenameMolecules(range, molname);
        _parser.IgnoreChilds();
    }
    if (el == "define") {
        string molname = attr["name"];
        string first = attr["first"];
        string nbeads = attr["nbeads"];
        string nmols = attr["nmols"];
        if (molname == "" && first == "" && nbeads == "" && nmols == "")
            throw runtime_error("invalid define tag");
        _top->CreateMoleculesByRange(molname,
                boost::lexical_cast<int>(first),
                boost::lexical_cast<int>(nbeads),
                boost::lexical_cast<int>(nmols));
        _parser.IgnoreChilds();
    }
}

void XMLTopologyReader::ParseBeadTypes(const string &el, map<string, string> &attr)
{
    cout << "XXX" << endl;
    if (el == "rename") {
        string name = attr["name"];
        string newname = attr["newname"];
        if (name == "" || newname == "")
            throw runtime_error("invalid rename tag");
        _top->RenameBeadType(name, newname);
        _parser.IgnoreChilds();
    }
}
}}

