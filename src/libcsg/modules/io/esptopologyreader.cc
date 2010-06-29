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
#include "esptopologyreader.h"

bool ESPTopologyReader::ReadTopology(string file, Topology &top)
{ 
    // cleanup topology to store new data
    ifstream fl;
    string parse_line, type, mass, tmp;
		int white_space1, white_space2, num_molecules, check_name;
		int *num_atoms = NULL;
		double box_x,box_y,box_z;
		bool HasMass;		
    top.Cleanup();

    fl.open(file.c_str());
    if(!fl.is_open())
        return false;
    
		// 'box_l' variable for box geometry
		////////////////////////////////////
		getline(fl, parse_line);
		parse_line = ReadBlockfileLine(parse_line, "box_l");
		white_space1 = parse_line.find(' ');
		if (white_space1 != string::npos) {
				white_space2 = parse_line.find(' ',white_space1+1);				
				if (white_space2 != string::npos) {
						box_x = atof(parse_line.substr(0,white_space1).c_str());
						box_y = atof(parse_line.substr(white_space1+1,white_space2).c_str());
						box_z = atof(parse_line.substr(white_space2).c_str());
				} else {
						cerr << "Error in parsing blockfile (box_l).\n";
						return false;						
				}
		} else {
				cerr << "Error in parsing blockfile (box_l).\n";
				return false;						
		}
		top.setBox(matrix(vec(box_x, 0, 0), vec(0, box_y, 0), vec(0, 0, box_z)));
		//cout << "Box geometry: " << box_x << " " << box_y << " " << box_z << endl;

		// 'num_molecules' tclvariable for total number of molecules
		////////////////////////////////////////////////////////////
    getline(fl, parse_line);
		parse_line = ReadBlockfileLine(parse_line, "num_molecules");
		num_molecules = atof(parse_line.c_str());
		//cout << "Number of molecules: " << num_molecules << endl;		

		// 'num_atoms' tclvariable for list of number of atoms (for each molecule)
		//////////////////////////////////////////////////////////////////////////
		num_atoms = (int*) calloc(num_molecules,sizeof(int));
		getline(fl, parse_line);
		parse_line = ReadBlockfileLine(parse_line, "num_atoms");
		for (int i=0; i<num_molecules-1;++i) {
				white_space1 = parse_line.find(' ');
				if (white_space1 != string::npos) {
						num_atoms[i] = atoi(parse_line.substr(0,white_space1+1).c_str());
						parse_line = parse_line.substr(white_space1+1);
				} else {
						cerr << "Error in parsing blockfile (num_atoms).\n";
						return false;
				}
		}
		num_atoms[num_molecules-1] = atoi(parse_line.c_str());		

		// 'particles' check whether mass is provided (if so, at the end)
		/////////////////////////////////////////////////////////////////
		getline(fl, parse_line);
		// make sure we're at 'particles'
		check_name = parse_line.find("{particles");
		if (check_name == string::npos) {
				cerr << "Can't find particles variable in blockfile.\n";
				return false;
		}
		check_name = parse_line.find("mass");
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
		
		// At this point we're assuming the particles have been ordered wrt molecules
		for (int mol=0; mol<num_molecules; ++mol) {
				Molecule *mi = top.CreateMolecule(boost::lexical_cast<string>(mol));
 				for (int atom=0; atom<num_atoms[mol]; ++atom) {
						fl >> tmp;	
						fl >> type;
						fl >> tmp;
						if (HasMass)
								fl >> mass;
						else
								mass = string("1.0");
						
						mi->AddBead(top.CreateBead(1, "", top.GetOrCreateBeadType(type), 0, atoi(mass.c_str()), 0),
												boost::lexical_cast<string>(atom));
 						getline(fl, parse_line);						
 				}
		}

		// Check that we've reached the end of the 'particles' variable
		getline(fl, parse_line);
		if (parse_line != "}") {
				cerr << "'num_molecules' and 'num_atoms' do not correspond to "
						"number of particles. Check .esp file.\n";
				return false;				
		}
		
		
		free(num_atoms);		
		
 		fl.close();
    
    cout << "WARNING: topology created from .esp file, charges aren't loaded!\n";
    
    return true;
}


string ESPTopologyReader::ReadBlockfileLine(string input_line, string variable)
{ 
		int pos_openvar, pos_closevar;		
		
		string str_find = "{" + variable;
		
		pos_openvar = input_line.find(str_find);
		if (pos_openvar == string::npos) {
				cerr << "Can't find '" << variable << "' variable in blockfile.\n";
				return false;
		}
		input_line = input_line.substr(pos_openvar+variable.length()+2);
		pos_closevar = input_line.find('}');
		if (pos_closevar == string::npos) {
				cerr << "Can't find '}' character in closing '" << variable << "' variable of the blockfile.\n";
				return false;
		}
		input_line = input_line.substr(0,pos_closevar);		
		
		return input_line;		
}
