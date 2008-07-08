#include "basis_set.h"
 /* chuncks of data which I will point the internal variables _nel_at and _nbasis_at and _basis_lbl_at too
   * this way I shall be able to deal with different basis sets
   */
void parce_string(string line, string delims, vector<string>*result);




static int      TZP_nel_at      [18] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
static int      TZP_nbasis_at   [18] = {6,6,20,20,20,20,20,20,20,20 ,27, 27, 27, 27, 27, 27, 27, 27 };  
static string   TZP_basis_lbl_1 [6]  = {"s","s","s","x","y","z"};
static string   TZP_basis_lbl_2 [20]  = {"s", "s", "s", "s", "s",
    					"x", "y", "z", "x", "y", "z", "x", "y", "z",	
					 "xx","yy", "zz", "xy", "xz", "yz"}; 
static string   TZP_basis_lbl_3 [27]  = {"s", "s", "s", "s", "s", "s",
    					"x","y","z", "x","y","z", "x","y","z", "x","y","z", "x","y","z",
					"xx","yy","zz","xy","xz","yz"};

static string*  TZP_basis_lbl_at[18] = {&TZP_basis_lbl_1[0],   
					&TZP_basis_lbl_1[0], 
					&TZP_basis_lbl_2[0], 
					&TZP_basis_lbl_2[0], 
					&TZP_basis_lbl_2[0], 
					&TZP_basis_lbl_2[0], 
					&TZP_basis_lbl_2[0], 
					&TZP_basis_lbl_2[0], 
					&TZP_basis_lbl_2[0], 
					&TZP_basis_lbl_2[0], 
					&TZP_basis_lbl_3[0], 
					&TZP_basis_lbl_3[0], 
					&TZP_basis_lbl_3[0], 
					&TZP_basis_lbl_3[0], 
					&TZP_basis_lbl_3[0], 
					&TZP_basis_lbl_3[0],
					&TZP_basis_lbl_3[0],
					&TZP_basis_lbl_3[0]};


static int      STO_nel_at      [18] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
static int      STO_nbasis_at   [18] = {1,1,5,5,5,5,5,5,5,5 ,9, 9, 9, 9, 9, 9, 9, 9 };
static string   STO_basis_lbl_1 [1]  = {"s"};
static string   STO_basis_lbl_2 [5]  = {"s", "s", "x", "y", "z"}; 
static string   STO_basis_lbl_3 [9]  = {"s", "s", "x", "y", "z","s" , "x", "y", "z"};
static string*  STO_basis_lbl_at[18] = {&STO_basis_lbl_1[0],   
					&STO_basis_lbl_1[0], 
					&STO_basis_lbl_2[0], 
					&STO_basis_lbl_2[0], 
					&STO_basis_lbl_2[0], 
					&STO_basis_lbl_2[0], 
					&STO_basis_lbl_2[0], 
					&STO_basis_lbl_2[0], 
					&STO_basis_lbl_2[0], 
					&STO_basis_lbl_2[0], 
					&STO_basis_lbl_3[0], 
					&STO_basis_lbl_3[0], 
					&STO_basis_lbl_3[0], 
					&STO_basis_lbl_3[0], 
					&STO_basis_lbl_3[0], 
					&STO_basis_lbl_3[0],
					&STO_basis_lbl_3[0],
					&STO_basis_lbl_3[0]};
					
static int      INDO_nel_at      [19] = {1,2,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,0};
static int      INDO_nbasis_at   [19] = {1,1,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,0};
static string   INDO_basis_lbl_1 [1]  = {"s"};
static string   INDO_basis_lbl_2 [4]  = {"s", "x", "y", "z"}; 
static string   INDO_basis_lbl_3 [4]  = {"s" , "x", "y", "z"};
static string*  INDO_basis_lbl_at[19] = {&INDO_basis_lbl_1[0],   
					&INDO_basis_lbl_1[0], 
					&INDO_basis_lbl_2[0], 
					&INDO_basis_lbl_2[0], 
					&INDO_basis_lbl_2[0], 
					&INDO_basis_lbl_2[0], 
					&INDO_basis_lbl_2[0], 
					&INDO_basis_lbl_2[0], 
					&INDO_basis_lbl_2[0], 
					&INDO_basis_lbl_2[0], 
					&INDO_basis_lbl_3[0], 
					&INDO_basis_lbl_3[0], 
					&INDO_basis_lbl_3[0], 
					&INDO_basis_lbl_3[0], 
					&INDO_basis_lbl_3[0], 
					&INDO_basis_lbl_3[0],
					&INDO_basis_lbl_3[0],
					&INDO_basis_lbl_3[0],
                                        &INDO_basis_lbl_1[0]};
					
void basis_set::set_basis_set(const string & a){
    if ( a == "INDO" ){
	_nel_at      = INDO_nel_at;
	_nbasis_at   = INDO_nbasis_at;
	_basis_lbl_at= INDO_basis_lbl_at;
    }
    else if ( a == "STO" ){
	_nel_at      = STO_nel_at;
	_nbasis_at   = STO_nbasis_at;
	_basis_lbl_at= STO_basis_lbl_at;
    }
    else if ( a == "TZP" ){
	_nel_at      = TZP_nel_at;
	_nbasis_at   = TZP_nbasis_at;
	_basis_lbl_at= TZP_basis_lbl_at;
    }
    else {
	cerr << a << " unknown basis set" <<endl;
    }
    basis_set_name = a;
}


int basis_set::set_bs_prim( const char * name_basis_file ){
    ifstream in(name_basis_file);
    if (!in){
	cerr << "File: " << name_basis_file << " non existant." <<endl;
	return 1;
    }
    string basis_name;
    in >> basis_name; 
    set_basis_set(basis_name);

    /* start reading in the info on primitives etc */
	
    string line;
    string delims=" \t";
    vector <string> words;
    int atom_now=0;
    while (in){
    	getline(in, line);
	words.clear();
	parce_string(line, delims, &words );
	if (words.size()==3){
		int lbl_at = atoi(words[1].c_str());
		int n_shells = atoi(words[2].c_str());
		if (lbl_at < atom_now ) {
		    cerr << "File: " << name_basis_file << " not in correct format " <<endl;
		    return  2;
		}
		while ( lbl_at > atom_now ){
			atom_shells empty;
			empty.n_shells =0;
			_atom_shells.push_back(empty);
			atom_now++;
		}
		int count_shells = 0;
		atom_shells full ;
		full.n_shells = n_shells;
		while ( true ) { // reading in the shells for atom whatever
			getline(in, line);
			words.clear();	
			parce_string(line, delims, &words );
			if (words.size() > 1 ) {   // readin in the specific shell for whatever
				count_shells = atoi(words[0].c_str()); // basically increment the counter (Im not chekcing)
				shell t_shell;
				t_shell._n_prim=0; // YES, its a  fuggin counter hence it start from uno
				if (words[1] == "s" ){
				    t_shell._lambda = 0;
				    t_shell._n_orbs = 1;
				}
				else if (words[1] == "p"){
				    t_shell._lambda = 1;
				    t_shell._n_orbs = 3;
				}
				else if (words[1] == "d"){
				    t_shell._lambda = 2;
				    t_shell._n_orbs = 6; // YES cartesian di functions: you got a problem wiht it punk?
				}
				else {
				    cerr << words[1] << " problem readin goddamned line: " << line << "wtf: order the ouput! (or expand the code)"<<endl;
				    return 3;
				}
				t_shell._n_prim++;
				t_shell._prim_exponent.push_back ( atof ( words[3].c_str() ) );
				t_shell._prim_contract.push_back ( atof ( words[4].c_str() ) );
				getline(in, line);
				words.clear();	
				parce_string(line, delims, &words );
				while ( words.size() > 1 ){
				    t_shell._n_prim++;
				    t_shell._prim_exponent.push_back ( atof ( words[3].c_str() ) );
				    t_shell._prim_contract.push_back ( atof ( words[4].c_str() ) );
				    getline(in, line);
				    words.clear();	
				    parce_string(line, delims, &words );
				}
				full._shells.push_back(t_shell);
				if (count_shells == n_shells ) break;
			}

		}
		_atom_shells.push_back(full);
		atom_now++;
	}
    }
    return 0;
}


void basis_set::print_all_primitive_info(ostream &out){
    for ( int j =0 ; j< _atom_shells.size(); ++j){
		
		for (int k=0; k< _atom_shells[j].n_shells;++k){
		    out << "Doing atom: " << j << " and shell: " << k <<endl;
		    out << _atom_shells[j]._shells[k]._n_orbs << endl;
		    out << _atom_shells[j]._shells[k]._lambda << endl;
		    out << _atom_shells[j]._shells[k]._n_prim << endl;
		    out << _atom_shells[j]._shells[k]._prim_exponent.size() <<endl;
		    out << _atom_shells[j]._shells[k]._prim_contract.size() <<endl;
		    for (int i=0; i<_atom_shells[j]._shells[k]._n_prim;++i){
			    out << _atom_shells[j]._shells[k]._prim_exponent[i] << endl;
			    out << _atom_shells[j]._shells[k]._prim_contract[i] << endl;
		    }
		}
	}
}
