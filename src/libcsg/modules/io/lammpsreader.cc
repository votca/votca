/* 
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#include <vector>
#include <boost/lexical_cast.hpp>
#include <votca/tools/getline.h>
#include "lammpsreader.h"

namespace votca { namespace csg {
    using namespace boost;
using namespace std;

bool LAMMPSReader::ReadTopology(string file,  Topology &top)
{
   _topology = true;
   top.Cleanup();

   _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open topology file: " + file);

   NextFrame(top);

    _fl.close();

    return true;
}

bool LAMMPSReader::Open(const string &file)
{
    _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open trajectory file: " + file);
    return true;
}

void LAMMPSReader::Close()
{
    _fl.close();
}

bool LAMMPSReader::FirstFrame(Topology &top)
{
    _topology = false;
    NextFrame(top);
    return true;
}

bool LAMMPSReader::NextFrame(Topology &top)
{
    string line;
    getline(_fl, line);
    while(!_fl.eof()) {
        if(line.substr(0, 5) != "ITEM:")
            throw std::ios_base::failure("unexpected line in lammps file:\n"+line);
        if(line.substr(6) == "TIMESTEP") {
                ReadTimestep(top, line);
        }
        else if(line.substr(6) == "NUMBER OF ATOMS") {
                ReadNumAtoms(top, line);
        }
        else if(line.substr(6) == "BOX BOUNDS") {
                ReadBox(top, line);
        }
        else if(line.substr(6, 5) == "ATOMS") {
                ReadAtoms(top, line);
                break;
        }

        else {
            throw std::ios_base::failure("unknown item lammps file : " + line.substr(6));
        }
        getline(_fl, line);
    }
    return !_fl.eof();;
}

void LAMMPSReader::ReadTimestep(Topology &top, string itemline)
{
    string s;
    getline(_fl, s);
    top.setStep(boost::lexical_cast<int>(s));
    cout << "Reading frame, timestep " << top.getStep() << endl;
}

void LAMMPSReader::ReadBox(Topology &top, string itemline)
{
    string s;

    matrix m;
    m.ZeroMatrix();
        
    for(int i=0; i<3; ++i) {
        getline(_fl, s);
        Tokenizer tok(s, " ");
        vector<double> v;
        tok.ConvertToVector(v);
        if(v.size()!=2)
            throw std::ios_base::failure("invalid box format");
        m[i][i] = v[1] - v[0];
    }
    top.setBox(m );
}

void LAMMPSReader::ReadNumAtoms(Topology &top, string itemline)
{
    string s;
    getline(_fl, s);
    _natoms = boost::lexical_cast<int>(s);
    if(!_topology && _natoms !=top.BeadCount())
        std::runtime_error("number of beads in topology and trajectory difffer");
}

void LAMMPSReader::ReadAtoms(Topology &top, string itemline) {
    top.CreateResidue("dum");
    bool pos=false;
    bool force=false;

    vector<string> fields;

    {
        Tokenizer tok(itemline.substr(12), " ");
        tok.ToVector(fields);
        for(Tokenizer::iterator i=tok.begin(); i!=tok.end(); ++i) {
            if(*i == "x" || *i == "y" || *i == "z")
                pos =true;
            else if(*i == "fx" || *i == "fy" || *i == "fz")
                force=true;
        }
    }
    for(int i=0; i<_natoms; ++i) {
        string s;
        getline(_fl, s);
        Bead *b;
        if(_topology)
            b = top.CreateBead(1, "", top.GetOrCreateBeadType("no"), 0, 0, 0);
        else
            b = top.getBead(i);

        b->HasPos(pos);
        b->HasF(force);

        Tokenizer tok(s, " ");
        Tokenizer::iterator itok= tok.begin();
        for(size_t j=0; itok!=tok.end(); ++itok, ++j) {
            if(j == fields.size())
                throw std::runtime_error("error, wrong number of columns in atoms section");
            else if(fields[j] == "x")
                b->Pos().x() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "y")
                b->Pos().y() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "z")
                b->Pos().z() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "fx")
                b->F().x() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "fy")
                b->F().y() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "fz")
                b->F().z() = boost::lexical_cast<double>(*itok);
            
        }
    }
}

}}

