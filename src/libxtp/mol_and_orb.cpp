/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/xtp/mol_and_orb.h>

namespace votca { namespace xtp {

map <string, int> MakePeriodicTable(){
    map <string, int> m;
    m.insert(pair<string, int>("H",0));
    m.insert(pair<string, int>("He",1));
    m.insert(pair<string, int>("Li",2));
    m.insert(pair<string, int>("Be",3));
    m.insert(pair<string, int>("B",4));
    m.insert(pair<string, int>("C",5));
    m.insert(pair<string, int>("N",6));
    m.insert(pair<string, int>("O",7));
    m.insert(pair<string, int>("F",8));
    m.insert(pair<string, int>("Ne",9));
    m.insert(pair<string, int>("Na",10));
    m.insert(pair<string, int>("Mg",11));
    m.insert(pair<string, int>("Al",12));
    m.insert(pair<string, int>("Si",13));
    m.insert(pair<string, int>("P",14));
    m.insert(pair<string, int>("S",15));
    m.insert(pair<string, int>("Cl",16));
    m.insert(pair<string, int>("Ar",17));
    m.insert(pair<string, int>("K",18));
    m.insert(pair<string, int>("Ca",19));
    m.insert(pair<string, int>("Sc",20));
    m.insert(pair<string, int>("Ti",21));
    m.insert(pair<string, int>("V",22));
    m.insert(pair<string, int>("Cr",23));
    m.insert(pair<string, int>("Mn",24));
    m.insert(pair<string, int>("Fe",25));
    m.insert(pair<string, int>("Co",26));
    m.insert(pair<string, int>("Ni",27));
    m.insert(pair<string, int>("Cu",28));
    m.insert(pair<string, int>("Zn",29));
    m.insert(pair<string, int>("Ga",30));
    m.insert(pair<string, int>("Ge",31));
    m.insert(pair<string, int>("As",32));
    m.insert(pair<string, int>("Se",33));
    m.insert(pair<string, int>("Br",34));
    m.insert(pair<string, int>("Kr",35));
    m.insert(pair<string, int>("Rb",36));
    m.insert(pair<string, int>("Sr",37));
    m.insert(pair<string, int>("Y",38));
    m.insert(pair<string, int>("Zr",39));
    m.insert(pair<string, int>("Nb",40));
    m.insert(pair<string, int>("Mo",41));
    m.insert(pair<string, int>("Tc",42));
    m.insert(pair<string, int>("Ru",43));
    m.insert(pair<string, int>("Rh",44));
    m.insert(pair<string, int>("Pd",45));
    m.insert(pair<string, int>("Ag",46));
    m.insert(pair<string, int>("Cd",47));
    m.insert(pair<string, int>("In",48));
    m.insert(pair<string, int>("Sn",49));
    m.insert(pair<string, int>("Sb",50));
    m.insert(pair<string, int>("Te",51));
    m.insert(pair<string, int>("I",52));
    m.insert(pair<string, int>("Xe",53));
    m.insert(pair<string, int>("Cs",54));
    m.insert(pair<string, int>("Ba",55));
    m.insert(pair<string, int>("La",56));
    m.insert(pair<string, int>("Ce",57));
    m.insert(pair<string, int>("Pr",58));
    m.insert(pair<string, int>("Nd",59));
    m.insert(pair<string, int>("Pm",60));
    m.insert(pair<string, int>("Sm",61));
    m.insert(pair<string, int>("Eu",62));
    m.insert(pair<string, int>("Gd",63));
    m.insert(pair<string, int>("Tb",64));
    m.insert(pair<string, int>("Dy",65));
    m.insert(pair<string, int>("Ho",66));
    m.insert(pair<string, int>("Er",67));
    m.insert(pair<string, int>("Tm",68));
    m.insert(pair<string, int>("Yb",69));
    m.insert(pair<string, int>("Lu",70));
    m.insert(pair<string, int>("Hf",71));
    m.insert(pair<string, int>("Ta",72));
    m.insert(pair<string, int>("W",73));
    m.insert(pair<string, int>("Re",74));
    m.insert(pair<string, int>("Os",75));
    m.insert(pair<string, int>("Ir",76));
    m.insert(pair<string, int>("Pt",77));
    m.insert(pair<string, int>("Au",78));
    m.insert(pair<string, int>("Hg",79));
    m.insert(pair<string, int>("Tl",80));
    m.insert(pair<string, int>("Pb",81));
    m.insert(pair<string, int>("Bi",82));
    m.insert(pair<string, int>("Po",83));
    m.insert(pair<string, int>("At",84));
    m.insert(pair<string, int>("Rn",85));
    return m ;
}

map <string, int> periodic_table = MakePeriodicTable();

/*void mol_and_orb::write_pdb(string file, string name_mol="PPY", const int & n=0 ){
    
    ofstream fl;
    static bool first_time=true;
    if (first_time==true){
       fl.open(file.c_str());
       first_time=false;
    }
    else {
        fl.open(file.c_str(), ios::app);
    }
    fl.setf(ios::fixed);
    for (int i =0 ; i < N ;i++ ){
         fl << setw(6) << "HETATM"
           << setw(5) << n+i
           << setw(5) << atom_labels[i]._type
           << setw(4) << name_mol
           << setw(3) << atom_labels[i]._lbl
           << setw(3) << "1" 
           << setw(4) << " "
           << setw(8) << setprecision(3) << RA*((atom_pos[i]).getX())
           << setw(8) << setprecision(3) << RA*((atom_pos[i]).getY())
           << setw(8) << setprecision(3) << RA*((atom_pos[i]).getZ())
           << setw(6) << setprecision(2) << 1.0
           << setw(6) << setprecision(2) << 20.0
           << setw(10) << "TM" << endl;
    }
    fl.close();
}*/

void mol_and_orb::write_pdb(string file, string name_mol, const int & n=0 ){
    FILE* out_qm;
    if(n==0){
       out_qm = fopen(file.c_str(), "w");
    }
    else{
        out_qm = fopen(file.c_str(), "a");
    }
    for (unsigned int i=0 ; i < N ;i++ ){
        fprintf(out_qm, "ATOM  %5d %4c %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", (n+i+1)%100000, (atom_labels[i]._type)[0],
name_mol.c_str(), " ", (n+1)/N, RA*((atom_pos[i]).getX()), RA*((atom_pos[i]).getY()), RA*((atom_pos[i]).getZ()), 1.0, 20.0);
    }
    fclose(out_qm);
}

int mol_and_orb::init_charges(const char * neutrfile, const char * crgfile){
    _neutr = new multipoles;
    _crged = new multipoles;
    _neutr->read_crg_eps(neutrfile);
    _crged->read_crg_eps(crgfile);
    if(_neutr->mpls.size() != this->N)
        throw std::runtime_error("number of atoms and entries in neutral charges file does not match");
    if(_crged->mpls.size() != this->N)
        throw std::runtime_error("number of atoms and entries in charged charges file does not match");
    return 0;
}

int mol_and_orb::init(const char * nameinput){
	int i=0, j=0 ;
	// read molecule//////////////////////////////////////////////////////////////////////////////////////

	ifstream in;
	in.open(nameinput);
	if(!in){cout << "Error, file " << nameinput << "does not exist" << endl; return 1;}
        
        // Skip header lines
        string line;
        std::getline(in, line);
        vector< string > split;
        Tokenizer toker(line, " \t");
        toker.ToVector(split);
        if (split.size() != 1) {
            cout << endl << "Error, file " << nameinput << " does not conform to xyz standards." << endl;
            throw runtime_error("(see above, input error)");
        }
        std::getline(in, line);
        
        // Continue with e x y z ...
	string  xch,ych,zch;
	double x,y,z;
        n_el=0;
	centre.setX(0.);
	centre.setY(0.);
	centre.setZ(0.);
        
	string word;
        
	vec pos;
	
	int lbl;
	string type;
        //char chtype;
	atom_type a_type;

	while ( in >> word){

		if(i%4==0){
		    if (word[0] != '#' ) {
                        type= word;
                        //chtype= word[0];
                    }
		    else {getline(in, word); i--;}
		}
              
		else if (i%4==1) xch = word;
		else if (i%4==2) ych = word;
		else if (i%4==3) {
			zch = word;
                        /*read the atom type*/
                        map<string, int>::iterator itp = periodic_table.find(type);
                        if (itp != periodic_table.end() ){
                            lbl = itp->second;
                            n_el += _bs->get_nel_at(lbl);
                            //cout << " Reading atom type " << type << ":" << itp->second << ":" << n_el << endl;
                        }
                        else {
                            cout << " Problem with atom type " << type << ":" << itp->second << ":" << n_el << endl;
                            throw runtime_error(string("Bad atom type, "));
                        }
			/*switch (chtype ){
			
			//set up labls, basis set ///////////////////////////////
			    case 'C': {  // definition of notation used:
				lbl = (5);
				n_el += _bs->get_nel_at(lbl);
				break;
				      }
			    case 'N': {
				lbl = ( 6);
				n_el += _bs->get_nel_at(lbl);
				break;
			
				      }
			    case 'O': {
				lbl = (7);
				n_el += _bs->get_nel_at(lbl);
				break;
				      }
			    case 'A': {  
				lbl = (12);
				n_el = _bs->get_nel_at(lbl);
				break;
				      }
			    case 'S': {
				lbl = (15);
				n_el += _bs->get_nel_at(lbl);
				break;
				      }
			    case 'H': {
				lbl = (0);
				n_el += _bs->get_nel_at(lbl);
				break; 
				      }
                            case 'B': {
                                lbl = (18);
				n_el += _bs->get_nel_at(lbl);
				break; 
                            }
			    default: {
				cout << "Bad atom type, we accept only: A (Al), C, O, H, N, S (use CAPITAL letters)" << endl;
				return 1;
				     }
			}*/
			////////////////////////////////////////////////////////
			sscanf(xch.c_str() , "%lf", &x );
			sscanf(ych.c_str() , "%lf", &y );
			sscanf(zch.c_str() , "%lf", &z );
			pos.setX(x / RA); pos.setY(y / RA) ; pos.setZ(z/RA);

			a_type._type = type;
			a_type._lbl  = lbl;

                        // cout << endl << "typ " << type << " lbl " << lbl << " " << x << " " << y << " " << z << endl;

			atom_labels.push_back(a_type);
			atom_pos.push_back(pos);
			
			j++;
		}
		i++;
	}
	N = j; 


	vector <vec> :: iterator iter_pos;
	for (iter_pos = atom_pos.begin() ;iter_pos < atom_pos.end() ; iter_pos++){
            centre += *iter_pos;
	}
	centre /= double(N);
        centre_input = centre;
	
//	for (iter_pos = atom_pos.begin() ;iter_pos < atom_pos.end() ; iter_pos++){
//	    *iter_pos -= centre;
//	}
//        centre.setX(0.); centre.setY(0.) ; centre.setZ(0.);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
        // cout << " Finished loading atoms" << endl;
	return 0;
}

int mol_and_orb::init_orbitals(orb & a, const char * namefile){
    vector <atom_type>:: iterator iter_at;
    int k=0;
    for (iter_at = atom_labels.begin() ; iter_at < atom_labels.end(); ++ iter_at){
    	k += _bs->get_nbasis_at( iter_at -> _lbl) ; 
    }
    string * basis;
    basis = new string[k];

    k=0;
    for (iter_at = atom_labels.begin() ; iter_at < atom_labels.end(); ++ iter_at){
    	for ( int i = 0 ; i < _bs->get_nbasis_at ( iter_at -> _lbl ) ;i++  ){
	   basis[k] = _bs->get_basis_lbl_at(iter_at -> _lbl , i) ;
	   k ++;
	}
    }
	
    a.init_orbitals( basis, k , namefile);
    orbitals = &a;
    delete [] basis;
    return 0;
}





void mol_and_orb::setAtomLabels(vector<string> types, vector<int> labels) {

    assert(types.size() == labels.size());

    for (unsigned int i = 0; i < types.size(); i++) {
        atom_type atom;
        atom._type = types[i];
        atom._lbl = labels[i];

        atom_labels.push_back(atom);


        map<string, int>::iterator itp = periodic_table.find(atom._type);
        if (itp != periodic_table.end() ){
            int lbl = itp->second;
            n_el += _bs->get_nel_at(lbl);
        }
        else {
            throw std::runtime_error("Unkown atom type " + atom._type);
        }
    }


    cout << "Number e-" << n_el << " --- Number atoms " << atom_labels.size() << endl;

}











}}
