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
        string xx = attr["xx"];
        string yy = attr["yy"];
        string zz = attr["zz"];
        if (xx == "" || yy == "" || zz == "")
            throw runtime_error("invalid box tag - missing xx, yy or zz");
        try {
           m[0][0] = boost::lexical_cast<double>(xx);
        } catch(boost::bad_lexical_cast &) {
           throw std::runtime_error("Cannot convert xx='"+ xx + "' to double");
        }
        try {
           m[1][1] = boost::lexical_cast<double>(yy);
        } catch(boost::bad_lexical_cast &) {
           throw std::runtime_error("Cannot convert yy='"+ yy + "' to double");
        }
        try {
           m[2][2] = boost::lexical_cast<double>(zz);
        } catch(boost::bad_lexical_cast &) {
           throw std::runtime_error("Cannot convert zz='"+ zz + "' to double");
        }
        _top->setBox(m);
        _parser.NextHandler(this, &XMLTopologyReader::ParseTopology);
    }
    else if(el == "beadtypes") {
         _parser.NextHandler(this, &XMLTopologyReader::ParseBeadTypes);
    } else if (el == "h5md_particle_group") {
        _top->setParticleGroup(attr["name"]);
        _parser.NextHandler(this, &XMLTopologyReader::ParseTopology);
    } else {
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
        if (molname == "" || first == "" || nbeads == "" || nmols == "")
            throw runtime_error("invalid define tag - missing molname, first, nbeads or nmols");
        int first_int, nbeads_int, nmols_int;
        try {
           first_int = boost::lexical_cast<int>(first);
        } catch(boost::bad_lexical_cast &) {
           throw std::runtime_error("Cannot convert first='"+ first + "' to int");
        }
        if (first_int < 1) {
           throw std::runtime_error("Attribute first is supose to be > 0, but found " + boost::lexical_cast<string>(first_int));
        }
        try {
           nbeads_int = boost::lexical_cast<int>(nbeads);
        } catch(boost::bad_lexical_cast &) {
           throw std::runtime_error("Cannot convert nbeads='"+ nbeads + "' to int");
        }
        if (nbeads_int < 1) {
           throw std::runtime_error("Attribute nbeads is supose to be > 0, but found " + boost::lexical_cast<string>(nbeads_int));
        }
        try {
           nmols_int = boost::lexical_cast<int>(nmols);
        } catch(boost::bad_lexical_cast &) {
           throw std::runtime_error("Cannot convert nmols='"+ nmols + "' to int");
        }
        if (nmols_int < 1) {
           throw std::runtime_error("Attribute nmols is supose to be > 0, but found " + boost::lexical_cast<string>(nmols_int));
        }
        _top->CreateMoleculesByRange(molname,first_int, nbeads_int, nmols_int);
        _parser.IgnoreChilds();
    }
}

void XMLTopologyReader::ParseBeadTypes(const string &el, map<string, string> &attr)
{
    if (el == "rename") {
        string name = attr["name"];
        string newname = attr["newname"];
        if (name == "" || newname == "")
            throw runtime_error("invalid rename tag");
        _top->RenameBeadType(name, newname);
        _parser.IgnoreChilds();
    }
    if (el == "mass") {
        string name = attr["name"];
        string svalue = attr["value"];
        if (name == "" || svalue == "")
            throw runtime_error("invalid mass tag");
	double value = boost::lexical_cast<double>(svalue);
        _top->SetBeadTypeMass(name, value);
        _parser.IgnoreChilds();
    }
}
}}

