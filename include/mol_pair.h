#ifndef MOL_PAIR
#define MOL_PAIR

#include <string>
#include "qm_molecule.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>

class mol_pair {

private:

	mol_and_orb *_morb1, *_morb2;

public:

	mol_pair(mol_and_orb* morb1, mol_and_orb* morb2) {
		_morb1 = morb1;
		_morb2 = morb2;
	}
	
	~mol_pair() {
		delete _morb1; // do we need this?
		delete _morb2;
	}	

	void write_conf(const char *file_type);
}

void mol_pair::write_conf(const char *file_type) {
	
	char atom_type;
	double x_coord, y_coord, z_coord;

	if ( !strcmp(file_type, "xyz") ) {
	//write xyz-file:
	// atom type: _morb1->atom_labels[i]._type;
	// atom x-coord: _morb1->atom_pos[i].x();
	// Num of atoms is 1 mol: _morb1->N;
		
		ofstream fout;
		fout.open("GEOM.xyz", ios::out);
		if (!fout) 
		    {
			cerr << "Opening file error!";
			exit(1);
		    }
		//writing num of atoms and a blank line
		fout << _morb1->N + _morb2->N << endl;
		fout << "\n" << endl;

		//writing first molecule:
		for (int i = 1; i <= _morb1->N; i++) {

		    atom_type =  _morb1->atom_labels[i]._type;
		    x_coord   =  _morb1->atom_pos[i].x();
		    y_coord   =  _morb1->atom_pos[i].y();
		    z_coord   =  _morb1->atom_pos[i].z();

		    fout << atom_type << "\t" << x_coord << "\t" << y_coord << "\t" << z_coord << endl; 
		}
		
		//writing second molecule:
		for (int i = 1; i <= _morb2->N; i++) {

		    atom_type =  _morb2->atom_labels[i]._type;
		    x_coord   =  _morb2->atom_pos[i].x();
		    y_coord   =  _morb2->atom_pos[i].y();
		    z_coord   =  _morb2->atom_pos[i].z();

		    fout << atom_type << "\t" << x_coord << "\t" << y_coord << "\t" << z_coord << endl; 		
		}
		
		fout.close();
	}
	
	else {
 		cout << "Unknown file format!" << endl;
	}

}

// int strcmp ( const char * str1, const char * str2 );
























#endif