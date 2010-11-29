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
#include "esptrajectoryreader.h"

namespace votca { namespace csg {

bool ESPTrajectoryReader::Open(const string &file)
{						
    _fl = file;		
    _temp_nextframe = string("tmp_nf_" + file);		
    return true;		
}


void ESPTrajectoryReader::Close()
{
    // Nothing to do - but define function since class herits from 
    // virtual class TrajectoryReader.
}


bool ESPTrajectoryReader::FirstFrame(Topology &top)
{
    ifstream file;
    ofstream out;		
    string parse_line, tmp;
    int atom;
    unsigned int check_name, check_name2;
    double r[3], v[3];		
    bool HasMass, HasVirtual;		

    file.open(_fl.c_str());

    if(!file.is_open())
	return false;		
		
    // Skip first 3 lines: variables 'box_l', 'num_molecules',
    // and 'num_atoms' read by ESPTopologyReader.
    getline(file, parse_line);
    getline(file, parse_line);
    getline(file, parse_line);

    // Check that the next line contains 'particles'
    getline(file, parse_line);
    check_name = parse_line.find("{particles");
    if (check_name == string::npos) {
	cerr << "Can't find particles variable in blockfile.\n";
	return false;
    }
    check_name  = parse_line.find("mass");
    check_name2 = parse_line.find("virtual");
    if (check_name != string::npos) {
	HasMass = 1;
	if (check_name2 != string::npos) {
	    HasVirtual = 1;
	    // make sure the format is correct
	    check_name = parse_line.find("{id type molecule mass virtual pos v}");
	    check_name2 = parse_line.find("{id type molecule mass virtual folded_position v}");
	} else {
	    HasVirtual = 0;
	    // make sure the format is correct
	    check_name = parse_line.find("{id type molecule mass pos v}");
	    check_name2 = parse_line.find("{id type molecule mass folded_position v}");
	}
    } else {
	HasMass = 0;
	if (check_name2 != string::npos) {
	    HasVirtual = 1;
	    // make sure the format is correct
	    check_name = parse_line.find("{id type molecule virtual pos v}");
	    check_name2 = parse_line.find("{id type molecule virtual folded_position v}");
	} else {
	    HasVirtual = 0;
	    // make sure the format is correct
	    check_name = parse_line.find("{id type molecule pos v}");
	    check_name2 = parse_line.find("{id type molecule folded_position v}");
	}
    }
    if (check_name == string::npos && check_name2 == string::npos) {
	cerr << "Check format of particles variable in blockfile.\n"
	    "Should be {id type molecule [mass] pos v}.\n"
	    "Instead: " << parse_line << endl;
	return false;
    }

    for (atom = 0; atom < top.BeadCount(); ++atom) {

	file >> tmp;
	file >> tmp;
	file >> tmp;				
	if (HasMass)
	    file >> tmp;				
	if (HasVirtual)
	    file >> tmp;

	// read particle position
	file >> r[0];
	file >> r[1];
	file >> r[2];
				
	// read particle velocities
	file >> v[0];
	file >> v[1];
	file >> v[2];				

	file >> tmp;				
				
	// Update particle properties
	top.getBead(atom)->setPos(r);
	top.getBead(atom)->setVel(v);

	getline(file,parse_line);
    }
    // Check that we've reached the end of the 'particles' variable
    getline(file, parse_line);
    if (parse_line != "}") {
	cerr << "'num_molecules' and 'num_atoms' do not correspond to "
	    "number of particles. Check .esp file.\n";
	return false;
    }

    getline(file, parse_line);		
    // Now that we've parsed the data for the first frame,
    // save the rest of the data into a new file (which will
    // be read by NextFrame(top).
    out.open(_temp_nextframe.c_str());	 
    while (!file.eof()) {
	out << parse_line << endl;
	getline(file, parse_line);
    }
		
    out.close();		
    file.close();		
    return true;
}

bool ESPTrajectoryReader::NextFrame(Topology &top)
{
    ifstream file;
    ofstream out;
    string parse_line, tmp;
    int atom; 
    unsigned int check_name;		
    bool HasMass;
    double r[3], v[3];		
		
		
    file.open(_temp_nextframe.c_str());

    if(!file.is_open())
        return false;
		
    if(file.eof()) return false;


    // Check that the next line contains 'particles'
    getline(file, parse_line);
    if (parse_line == "") return false;		
    check_name = parse_line.find("{particles");
    if (check_name == string::npos) {
	cerr << "Can't find particles variable in blockfile.\n";
	return false;
    }
    if (check_name != string::npos) {
	HasMass = 1;
	// make sure the format is correct
	check_name = parse_line.find("{id type molecule mass pos v}");
    } else {
	HasMass = 0;
	// make sure the format is correct
	check_name = parse_line.find("{id type molecule pos v}");
    }
    if (check_name == string::npos) {
	cerr << "Check format of particles variable in blockfile.\n"
	    "Should be {id type molecule [mass] pos v}.\n"
	    "Instead: " << parse_line << endl;
	return false;
    }

    for (atom = 0; atom < top.BeadCount(); ++atom) {

	file >> tmp;
	file >> tmp;
	file >> tmp;				
	if (HasMass)
	    file >> tmp;				

	// read particle position
	file >> r[0];
	file >> r[1];
	file >> r[2];
				
	// read particle velocities
	file >> v[0];
	file >> v[1];
	file >> v[2];				

	file >> tmp;				
				
	// Update particle properties
	top.getBead(atom)->setPos(r);
	top.getBead(atom)->setVel(v);

	getline(file,parse_line);
    }
    // Check that we've reached the end of the 'particles' variable
    getline(file, parse_line);
    if (parse_line != "}") {
	cerr << "'num_molecules' and 'num_atoms' do not correspond to "
	    "number of particles. Check .esp file.\n";
	return false;
    }

    getline(file, parse_line);		
    // Now that we've parsed the data for the first frame,
    // save the rest of the data into a new file (which will
    // be read by NextFrame(top).
    out.open(_temp_nextframe.c_str());	 
    while (!file.eof()) {
	out << parse_line << endl;
	getline(file, parse_line);
    }
		
    out.close();		
    file.close();		

    return true;
}

}}

