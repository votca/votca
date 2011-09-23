/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#ifndef BASIS_SET_H
#define BASIS_SET_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "votca/tools/tokenizer.h"
#include <boost/lexical_cast.hpp>
#include <votca/tools/property.h>

namespace votca { namespace moo {

using namespace std;
using namespace boost;
using namespace votca::tools;

struct shell{
	    int _n_orbs;
	    int _lambda;
	    int _n_prim;
	    vector <double> _prim_exponent;
	    vector <double> _prim_contract;
            shell() : _n_orbs(0), _lambda(0), _n_prim(0) {}
	};

class basis_set{
    private:
	
	struct atom_shells{
	    int n_shells;
	    vector <shell> _shells;
	};
	int    *  _nel_at;
	int    *  _nbasis_at;
	string ** _basis_lbl_at;

	string basis_set_name;
        Property _options; // to parse xml files

	vector < atom_shells >   _atom_shells;

        void parse_xml_basisset_info(const string &);
        //void ParseAtomBasisSet(xmlDocPtr doc, xmlNodePtr cur);

    public:
	// defaul constructor defaults to INDO
	basis_set() {
		set_basis_set("INDO");
	}
	~basis_set(){
            
        #ifdef DEBUG
        cout << "callgin the nasis set classdestructor" << endl;
        #endif
        }

	basis_set & operator =(const basis_set a);

	const int & get_nel_at(const int & i) const {
		return _nel_at[i];
	}

	const int & get_nbasis_at(const int & i) const {
		return _nbasis_at[i];
	} 

	const string & get_basis_lbl_at(const int & i, const int & j) const {
		return _basis_lbl_at[i][j];
	} 

	const string & get_basis_set_name() const {
	    return basis_set_name;
	}

        void set_basis_set(const string & a); 
	int set_bs_prim(const char * ) ;

	/// return the number of shells for atom of lbl i
	const int & get_n_shells( const int & i ) const {
		return _atom_shells[i].n_shells;
	}
	
	const shell * get_shell(const int &i, const int &j ) const {
		return &(_atom_shells[i]._shells[j]);
	}
	 

	/// return the number of primitives in the jth shell of the ith atom
	const int & get_n_primitives( const int & i, const int & j) const {
		return _atom_shells[i]._shells[j]._n_prim;
	}

	/// returns the number of basis sets of the jth shell on the ith atom
	const int & get_n_basis_shell(const int & i, const int & j ) const {
	    return _atom_shells[i]._shells[j]._n_orbs;
	}

	///return the angular momentum of the jth shell on the ith atom
	const int & get_am_shell(const int &i, const int & j) const {
		return _atom_shells[i]._shells[j]._lambda;
	}

	/// return the exponent of the kth primitve on the jth shell on the ith atom
	const double &get_exponent_prim(const int &i, const int &j, const int &k ){
	    return _atom_shells[i]._shells[j]._prim_exponent[k];
	}
	
	/// return the exponent of the kth primitve on the jth shell on the ith atom
	const double &get_contract_prim(const int &i, const int &j, const int &k ){
	    return _atom_shells[i]._shells[j]._prim_contract[k];
	}

	void print_all_primitive_info(ostream &out);
};

}}

#endif
