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

#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <votca/tools/getline.h>
#include <boost/algorithm/string.hpp>
#include "groreader.h"

namespace votca { namespace csg {

bool GROReader::ReadTopology(string file,  Topology &top)
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

bool GROReader::Open(const string &file)
{
    _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open trajectory file: " + file);
    return true;
}

void GROReader::Close()
{
    _fl.close();
}

bool GROReader::FirstFrame(Topology &top)
{
    _topology = false;
    NextFrame(top);
    return true;
}

bool GROReader::NextFrame(Topology &top)
{
    string tmp;
    getline(_fl, tmp);//title
    if (_fl.eof()) {
      return !_fl.eof();
    }
    getline(_fl, tmp);//number atoms
    int natoms = atoi(tmp.c_str());
    if(!_topology && natoms != top.BeadCount())
      throw std::runtime_error("number of beads in topology and trajectory differ");
   
    for(int i=0;i<natoms; i++) {
        string line;
        getline(_fl, line);
	string resNum,resName, atName, x,y,z;
	try {
          resNum= string(line,0,5); // %5i
          resName = string(line,5,5); //%5s
	  atName = string(line,10,5); // %5s
	  //atNum= string(line,15,5); // %5i not needed
	  x = string(line,20,8); // %8.3f
	  y = string(line,28,8); // %8.3f
	  z = string(line,36,8); // %8.3f
	} catch (std::out_of_range& err) {
	  throw std::runtime_error("Misformated gro file");
	}
        boost::algorithm::trim(atName);
        boost::algorithm::trim(resName);
        boost::algorithm::trim(resNum);
        boost::algorithm::trim(x);
        boost::algorithm::trim(y);
        boost::algorithm::trim(z);
	string vx,vy,vz;
	bool hasVel=true;
	try {
	  vx = string(line,44,8); // %8.4f
	  vy = string(line,52,8); // %8.4f
	  vz = string(line,60,8); // %8.4f
	} catch (std::out_of_range& err) {
	  hasVel=false;
	}

        Bead *b;
        if(_topology){
	  int resnr = boost::lexical_cast<int>(resNum);
          if (resnr < 1)
            throw std::runtime_error("Misformated gro file, resnr has to be > 0");
	  //TODO: fix the case that resnr is not in ascending order
	  if(resnr > top.ResidueCount()) {
	    while ((resnr-1)>top.ResidueCount()){ //gro resnr should start with 1 but accept sloppy files
	      top.CreateResidue("DUMMY"); // create dummy residue, hopefully it will never show
	      cout << "Warning: residue numbers not continous, create DUMMY residue with nr " << top.ResidueCount() << endl;
	    }
            top.CreateResidue(resName);
	  }
          //this is not correct, but still better than no type at all!
	  auto type = top.GetOrCreateBeadType(atName);

	  // res -1 as internal number starts with 0
	  b = top.CreateBead(1, atName, &type, resnr-1, 1., 0.);
	} else {
          b = top.getBead(i);
	}

        b->setPos(vec(
          boost::lexical_cast<double>(x),
          boost::lexical_cast<double>(y),
          boost::lexical_cast<double>(z)
        ));
	if (hasVel) {
          boost::algorithm::trim(vx);
          boost::algorithm::trim(vy);
          boost::algorithm::trim(vz);
          b->setVel(vec(
            boost::lexical_cast<double>(vx),
            boost::lexical_cast<double>(vy),
            boost::lexical_cast<double>(vz)
          ));
	}
    }

    getline(_fl, tmp); //read box line
    if(_fl.eof())
      throw std::runtime_error("unexpected end of file in poly file, when boxline");
    Tokenizer tok(tmp, " ");
    vector<double> fields;
    tok.ConvertToVector<double>(fields);
    matrix box;
    if ( fields.size() == 3 ) {
      box.ZeroMatrix();
      for (int i=0;i<3;i++) box[i][i] = fields[i];
    } else if ( fields.size() == 9 ) {
      box[0][0]=  fields[0];
      box[1][1]=  fields[1];
      box[2][2]=  fields[2];
      box[1][0]=  fields[3];
      box[2][0]=  fields[4];
      box[0][1]=  fields[5];
      box[2][1]=  fields[6];
      box[0][2]=  fields[7];
      box[1][2]=  fields[8];
    } else {
      throw std::runtime_error("Error while reading box (last) line");
    }
    top.setBox(box);

    if (_topology) {
      cout << "WARNING: topology created from .gro file, masses and charges are wrong!\n";
    }
    
    return !_fl.eof();
}

}}
