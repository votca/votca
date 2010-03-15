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

#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <stdio.h>
#include <stdexcept>
#include <expat.h>
#include "xmltopologyreader.h"

static void start_hndl(void *data, const char *el, const char **attr)
{
    XMLTopologyReader *reader =
        (XMLTopologyReader*)XML_GetUserData((XML_Parser*)data);

    map<string, string> mattr;

    for (int i = 0; attr[i]; i += 2)
       mattr[attr[i]] = attr[i + 1];
    string sel = el;
    reader->StartElemHndl(sel, mattr);
}

static void end_hndl(void *data, const char *el)
{
    XMLTopologyReader *reader =
        (XMLTopologyReader*)XML_GetUserData((XML_Parser*)data);
    reader->EndElemHndl(el);
}

bool XMLTopologyReader::ReadTopology(string filename, Topology &top)
{ 
  _top = &top;
  
  XML_Parser parser = XML_ParserCreate(NULL);
  if (! parser)
    throw std::runtime_error("Couldn't allocate memory for xml parser");

  XML_UseParserAsHandlerArg(parser);
  XML_SetElementHandler(parser, start_hndl, end_hndl);
//  XML_SetCharacterDataHandler(parser, char_hndl);

  ifstream fl;
  fl.open(filename.c_str());
  if(!fl.is_open())
    throw std::ios_base::failure("Error on open xml topology: " + filename);
  
  SetHandler(&XMLTopologyReader::ParseRoot);
  XML_SetUserData(parser, (void*)this);
  while(!fl.eof()) {
    string line;
    getline(fl, line);
    line=line + "\n";
    if (! XML_Parse(parser, line.c_str(), line.length(), false))
      throw  std::ios_base::failure(filename + ": Parse error at line " +
          boost::lexical_cast<string>(XML_GetCurrentLineNumber(parser)) + "\n" +
          XML_ErrorString(XML_GetErrorCode(parser)));
  }
  fl.close();
    
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
        SetHandler(&XMLTopologyReader::ParseTopology);
    }
    else {
        throw std::runtime_error("wrong root node in xml topology file");
    }
}

void XMLTopologyReader::ParseTopology(const string &el, map<string, string> &attr)
{
    if(el == "molecules")
        SetHandler(&XMLTopologyReader::ParseMolecules);
}

void XMLTopologyReader::ParseMolecules(const string &el, map<string, string> &attr)
{
    map<string, string>::iterator iter;
    if (el == "clear") {
        _top->ClearMoleculeList();
        SetHandler(&XMLTopologyReader::ParseIgnore);
    }
    else if (el == "rename") {
        string molname = attr["name"];
        string range = attr["range"];
        if (molname == "" || range == "")
            throw runtime_error("invalid rename tag");
        _top->RenameMolecules(range, molname);
        SetHandler(&XMLTopologyReader::ParseIgnore);
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
        SetHandler(&XMLTopologyReader::ParseIgnore);
    }
}

void XMLTopologyReader::ParseIgnore(const string &el, map<string, string> &attr)
{
    SetHandler(&XMLTopologyReader::ParseIgnore);
}

void XMLTopologyReader::StartElemHndl(const string &el, map<string, string> &attr)
{
    (this->*_handler)(el, attr);
}

void XMLTopologyReader::EndElemHndl(const string &el)
{
    _stack_handler.pop();
    _handler = _stack_handler.top();
}

void XMLTopologyReader::SetHandler(ElemHandler_t handler)
{
    _handler = handler;
    _stack_handler.push(handler);
}
