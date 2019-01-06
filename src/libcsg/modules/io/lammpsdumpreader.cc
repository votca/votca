/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
#include "lammpsdumpreader.h"

namespace votca { namespace csg {
    using namespace boost;
using namespace std;

bool LAMMPSDumpReader::ReadTopology(string file,  Topology &top)
{
   _topology = true;
   top.Cleanup();

   _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open topology file: " + file);
   _fname=file;

   NextFrame(top);

    _fl.close();

    return true;
}

bool LAMMPSDumpReader::Open(const string &file)
{
    _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open trajectory file: " + file);
    _fname=file;
    return true;
}

void LAMMPSDumpReader::Close()
{
    _fl.close();
}

bool LAMMPSDumpReader::FirstFrame(Topology &top)
{
    _topology = false;
    NextFrame(top);
    return true;
}

bool LAMMPSDumpReader::NextFrame(Topology &top)
{
    string line;
    getline(_fl, line);
    while(!_fl.eof()) {
        if(line.substr(0, 5) != "ITEM:")
            throw std::ios_base::failure("unexpected line in lammps file:\n"+line);
        if(line.substr(6,8) == "TIMESTEP") {
                ReadTimestep(top, line);
        }
        else if(line.substr(6,15) == "NUMBER OF ATOMS") {
                ReadNumAtoms(top, line);
        }
        else if(line.substr(6,10) == "BOX BOUNDS") {
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
    if (_topology) {
      cout << "WARNING: topology created from .dump file, masses, charges, types, residue names are wrong!\n";
    }
    return !_fl.eof();;
}

void LAMMPSDumpReader::ReadTimestep(Topology &top, string itemline)
{
    string s;
    getline(_fl, s);
    top.setStep(boost::lexical_cast<int>(s));
    cout << "Reading frame, timestep " << top.getStep() << endl;
}

void LAMMPSDumpReader::ReadBox(Topology &top, string itemline)
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

void LAMMPSDumpReader::ReadNumAtoms(Topology &top, string itemline)
{
    string s;
    getline(_fl, s);
    _natoms = boost::lexical_cast<int>(s);
    if(!_topology && _natoms !=top.BeadCount())
        std::runtime_error("number of beads in topology and trajectory differ");
}

void LAMMPSDumpReader::ReadAtoms(Topology &top, string itemline) {
    if(_topology){
      top.CreateResidue("dum");
      for(int i=0; i<_natoms; ++i) {
        (void)top.CreateBead(1, "no", (top.GetOrCreateBeadType("no")), 0, 0, 0);
      }
    }

    bool pos=false;
    bool force=false;
    bool vel=false;
    int id=-1;

    vector<string> fields;

    {
        Tokenizer tok(itemline.substr(12), " ");
        tok.ToVector(fields);
	int j=0;
        for(Tokenizer::iterator i=tok.begin(); i!=tok.end(); ++i,++j) {
            if(*i == "x" || *i == "y" || *i == "z")
                pos =true;
	    else if(*i == "xu" || *i == "yu" || *i == "zu")
                pos =true;
	    else if(*i == "xs" || *i == "ys" || *i == "zs")
                pos =true;
	    else if(*i == "vx" || *i == "vy" || *i == "vz")
                vel =true;
            else if(*i == "fx" || *i == "fy" || *i == "fz")
                force=true;
	    else if(*i == "id")
	      id=j;
        }
    }
    if (id<0){
      throw std::runtime_error("error, id not found in any column of the atoms section");
    }

    for(int i=0; i<_natoms; ++i) {
        string s;
        getline(_fl, s);
        if (_fl.eof())
           throw std::runtime_error("Error: unexpected end of lammps file '" + _fname + "' only " + boost::lexical_cast<string>(i) + " atoms of " + boost::lexical_cast<string>(_natoms) + " read.");

        Tokenizer tok(s, " ");
        Tokenizer::iterator itok= tok.begin();
        vector<string> fields2;
        tok.ToVector(fields2);
	// internal numbering begins with 0
	int atom_id = boost::lexical_cast<int>(fields2[id]);
	if (atom_id > _natoms)
           throw std::runtime_error("Error: found atom with id "+ boost::lexical_cast<string>(atom_id) + " but only "+ boost::lexical_cast<string>(_natoms) + " atoms defined in header of file '" + _fname + "'");
        Bead *b = top.getBead(atom_id-1);
        b->HasPos(pos);
        b->HasF(force);
        b->HasVel(vel);
        matrix m=top.getBox();

        for(size_t j=0; itok!=tok.end(); ++itok, ++j) {
            if(j == fields.size())
                throw std::runtime_error("error, wrong number of columns in atoms section");
            else if(fields[j] == "x")
                b->Pos().x() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "y")
                b->Pos().y() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "z")
                b->Pos().z() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "xu")
                b->Pos().x() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "yu")
                b->Pos().y() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "zu")
                b->Pos().z() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "xs")
                b->Pos().x() = boost::lexical_cast<double>(*itok)*m[0][0];
            else if(fields[j] == "ys")
                b->Pos().y() = boost::lexical_cast<double>(*itok)*m[1][1];
            else if(fields[j] == "zs")
                b->Pos().z() = boost::lexical_cast<double>(*itok)*m[2][2];
            else if(fields[j] == "vx")
                b->Vel().x() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "vy")
                b->Vel().y() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "vz")
                b->Vel().z() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "fx")
                b->F().x() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "fy")
                b->F().y() = boost::lexical_cast<double>(*itok);
            else if(fields[j] == "fz")
                b->F().z() = boost::lexical_cast<double>(*itok);
            else if((fields[j] == "type")&&_topology){
                auto type = top.GetOrCreateBeadType(*itok);
                b->setType(type);
            }
        }
    }
}

}}

