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

#include <cstdlib>
#include <iostream>
#include <votca/csg/topology.h>
#include "dlpolytrajectoryreader.h"

namespace votca { namespace csg {

using namespace std;

bool DLPOLYTrajectoryReader::Open(const string &file)
{
    _fl.open("HISTORY");
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open topology file: HISTORY");
    return true;
}

void DLPOLYTrajectoryReader::Close()
{
    _fl.close();
}

bool DLPOLYTrajectoryReader::FirstFrame(Topology &top)
{
    _first_frame=true;
    bool res=NextFrame(top);
    _first_frame=false;
    return res;
}

bool DLPOLYTrajectoryReader::NextFrame(Topology &top)
{
    string line;
    if (_first_frame) {
      getline(_fl, line); //title
      getline(_fl, line); // 2nd header line
      Tokenizer tok(line, " ");
      vector<int> fields;
      tok.ConvertToVector<int>(fields);
      if(fields[2] !=top.BeadCount())
        throw std::runtime_error("number of beads in topology and trajectory differ");
    }

    //read normal frame
    getline(_fl, line); // timestep line
    if(!_fl.eof()) {
      int natoms;
      {
        Tokenizer tok(line, " ");
        vector<string> fields;
        tok.ToVector(fields);
        natoms = boost::lexical_cast<int>(fields[2]);
        if(natoms !=top.BeadCount())
          throw std::runtime_error("number of beads in topology and trajectory differ");
	top.setTime(boost::lexical_cast<double>(fields[6]));
      }

      vec box_vectors[3];
      for (int i=0;i<3;i++){ // read 3 box lines
        getline(_fl, line);
        if(_fl.eof())
          throw std::runtime_error("unexpected end of file in dlpoly file, when reading box vector"  +
	      boost::lexical_cast<string>(i));
        Tokenizer tok(line, " ");
        vector<double> fields;
        tok.ConvertToVector<double>(fields);
	box_vectors[i]=vec(fields[0],fields[1],fields[2]);
      }
      matrix box(box_vectors[0],box_vectors[1],box_vectors[2]);
      top.setBox(box);

      for (int i=0;i<natoms;i++){

	{
          getline(_fl, line); //atom header line
          if(_fl.eof())
            throw std::runtime_error("unexpected end of file in dlpoly filei when reader atom nr" +
		boost::lexical_cast<string>(i+1));

          vector<string> fields;
          Tokenizer tok(line, " ");
          tok.ToVector(fields);
	  int id=boost::lexical_cast<double>(fields[1]);
	  if (i+1 != id )
            throw std::runtime_error("unexpected atom number in dlpoly file, expected" +
		 boost::lexical_cast<string>(i+1) + "got" + boost::lexical_cast<string>(id));
	}

	Bead *b = top.getBead(i);
	vec atom_vec[3];
	for (int j=0;j<3;j++){
          getline(_fl, line); //read atom positions
          if(_fl.eof())
            throw std::runtime_error("unexpected end of file in dlpoly file when reading atom vector" +
		boost::lexical_cast<string>(j) + " of atom " + boost::lexical_cast<string>(i+1) );
          vector<double> fields;
          Tokenizer tok(line, " ");
          tok.ConvertToVector<double>(fields);
	  atom_vec[j]=vec(fields[0],fields[1],fields[2]);
	}
	b->setPos(atom_vec[0]);
	b->setVel(atom_vec[1]);
	b->setF(atom_vec[2]);
      }
    }
    return !_fl.eof();
}

}}
