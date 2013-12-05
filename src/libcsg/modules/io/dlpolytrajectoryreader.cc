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
#include <boost/filesystem/convenience.hpp> 
#include <votca/csg/topology.h>
#include <votca/csg/boundarycondition.h>
#include <votca/tools/getline.h>
#include "dlpolytrajectoryreader.h"

namespace votca { namespace csg {

using namespace std;

bool DLPOLYTrajectoryReader::Open(const string &file)
{
    boost::filesystem::path filepath(file.c_str());
    string filename;

    if ( boost::filesystem::basename(filepath).size() == 0 ) {
      if (filepath.parent_path().string().size() == 0) {
        filename="HISTORY";
      } else {
	filename=filepath.parent_path().string() + "/HISTORY";
      }
    } else {
      filename=file;
    }

    _fl.open(filename.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open trajectory file: "+ filename);
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
    static int  mavecs=0;
    static int  mpbct=0;
    static int  matoms=0;
    static bool hasVs=false;
    static bool hasFs=false;

    string line;

    BoundaryCondition::eBoxtype pbc_type=BoundaryCondition::typeAuto;

    if (_first_frame) {
      getline(_fl, line); //title
      cout << "Read from HISTORY: '" << line << "' - header" << endl;

      getline(_fl, line); // 2nd header line
      cout << "Read from HISTORY: '" << line << "' - directives line" << endl;

      Tokenizer tok(line, " \t");
      vector<int> fields;
      tok.ConvertToVector<int>(fields);

      mavecs = boost::lexical_cast<int>(fields[0]);
      mpbct  = boost::lexical_cast<int>(fields[1]);
      matoms = boost::lexical_cast<int>(fields[2]);

      //if(mavecs > 0 ) hasVs=true;
      //if(mavecs > 1 ) hasFs=true;
      hasVs=(mavecs > 0);
      hasFs=(mavecs > 1);

      if(hasVs != top.HasVel() || hasFs != top.HasForce()) {
	//cout << "Warning: TrajKey (# of atom vectors) in HISTORY (trajectory) header & CONFIG (initial frame) differ" << endl;
	top.SetHasVel(hasVs);
	top.SetHasForce(hasFs);
      }

      //if(fields[2] != top.BeadCount())
      if(matoms != top.BeadCount())
        throw std::runtime_error("Number of atoms/beads in HISTORY (trajectory) header & CONFIG (initial frame) differ");

      if(mpbct == 0) {
	pbc_type=BoundaryCondition::typeOpen;
      }
      else if(mpbct == 1 || mpbct == 2 ) {
	pbc_type=BoundaryCondition::typeOrthorhombic;
      }
      else if(mpbct == 3) {
	pbc_type=BoundaryCondition::typeTriclinic;
      }

      if(pbc_type != top.getBoxType())
	//throw std::runtime_error("PBC type in HISTORY (trajectory) header & CONFIG (initial frame) differ");
	cout << "Warning: PBC type in HISTORY (trajectory) header & CONFIG (initial frame) differ" << endl;
    }

    //read normal frame
    getline(_fl, line); // timestep line
    cout << "Read from HISTORY: '" << line << "'" << endl;

    if(!_fl.eof()) {
      double dtime,stime,scale;
      int nstep;
      int natoms;
      int navecs;
      int npbct;

      //scale = 0.1; // AB: factor to convert Angstroem to nm
      scale = 1.0; // AB: factor to convert Angstroem to nm - not needed?

      {
        Tokenizer tok(line, " \t");
        vector<string> fields;
        tok.ToVector(fields);

	nstep  = boost::lexical_cast<int>(fields[1]);
        natoms = boost::lexical_cast<int>(fields[2]);
	navecs = boost::lexical_cast<int>(fields[3]);
	npbct  = boost::lexical_cast<int>(fields[4]);
	dtime  = boost::lexical_cast<double>(fields[5]);
	stime  = boost::lexical_cast<double>(fields[fields.size()-1]);

	//cout << "Read from HISTORY: natoms = " << natoms << ", levcfg = " << fields[3];
	//cout << ", dt = " << fields[5] << ", time = " << stime  << endl;

	// AB-TODO: would be nice to store navecs (levcfg) & npbct (imcon, PBC switch)

        if(natoms != top.BeadCount())
          throw std::runtime_error("number of atoms/beads in HISTORY (trajectory) frame & FIELD (topology) differ");
        if(navecs != mavecs)
	  throw std::runtime_error("TrajKey (# of atom vectors) in HISTORY (trajectory) header & frame differ");
        if(natoms != matoms)
	  throw std::runtime_error("Number of atoms/beads in HISTORY (trajectory) header & frame differ");
        if(npbct != mpbct)
	  throw std::runtime_error("PBC type in HISTORY (trajectory) header & frame differ");

	top.setTime(nstep*dtime);
	//top.setTime(boost::lexical_cast<double>(fields[5]));

	if(stime != top.getTime() )
	  cout << "Check: nstep = " << nstep << ", dt = " << dtime << ", time = " << top.getTime() << " (correct?)" << endl;
      }

      vec box_vectors[3];
      for (int i=0;i<3;i++){ // read 3 box lines
        getline(_fl, line);
	//cout << "Read from HISTORY: '" << line << "' - box vector # " << i+1 << endl;

        if(_fl.eof())
          throw std::runtime_error("unexpected EOF in HISTOPY (trajectory), when reading box vector"+
	      boost::lexical_cast<string>(i));
	
        Tokenizer tok(line, " \t");
        vector<double> fields;
        tok.ConvertToVector<double>(fields);
	box_vectors[i]=vec(fields[0]*scale,fields[1]*scale,fields[2]*scale);
      }
      matrix box(box_vectors[0],box_vectors[1],box_vectors[2]);

      if(npbct == 0) {
	pbc_type=BoundaryCondition::typeOpen;
      }
      else if(npbct == 1 || npbct == 2 ) {
	pbc_type=BoundaryCondition::typeOrthorhombic;
      }
      else if(npbct == 3) {
	pbc_type=BoundaryCondition::typeTriclinic;
      }

      top.setBox(box,pbc_type);

      for (int i=0;i<natoms;i++){

	{
          getline(_fl, line); //atom header line
	  //cout << "Read from HISTORY: '" << line << endl;

          if(_fl.eof())
            throw std::runtime_error("unexpected EOF in HISTOPY (trajectory), when reading atom/bead nr" +
		boost::lexical_cast<string>(i+1));

          vector<string> fields;
          Tokenizer tok(line, " \t");
          tok.ToVector(fields);
	  int id=boost::lexical_cast<double>(fields[1]);
	  if (i+1 != id )
            throw std::runtime_error("unexpected atom/bead index in HISTOPY (trajectory), expected" +
		 boost::lexical_cast<string>(i+1) + "got" + boost::lexical_cast<string>(id));
	}

	Bead *b = top.getBead(i);
	vec atom_vec[3];
	for (int j=0;j<min(navecs,2)+1;j++){
          getline(_fl, line); //read atom positions
          if(_fl.eof())
            throw std::runtime_error("unexpected EOF in HISTOPY (trajectory) when reading atom/bead vector" +
		boost::lexical_cast<string>(j) + " of atom " + boost::lexical_cast<string>(i+1) );
          vector<double> fields;
          Tokenizer tok(line, " \t");
          tok.ConvertToVector<double>(fields);
	  if(j<3) { 
	    // AB: convert Angstroem to nm; j=1 => x,y,z; j=2 => vx,vy,vz
	    atom_vec[j]=vec(fields[0]*scale,fields[1]*scale,fields[2]*scale);
	  }
	  else { // AB: forces do not need to be converted (assuming kJ/mole?)
	    atom_vec[j]=vec(fields[0],fields[1],fields[2]);
	  }
	}
	b->setPos(atom_vec[0]);
	if (navecs > 0)
	  b->setVel(atom_vec[1]);
	if (navecs > 1)
	  b->setF(atom_vec[2]);
      }
    }
    return !_fl.eof();
}

}}
