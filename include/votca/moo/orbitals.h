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

#ifndef FILE_ORBITALS_H
#define FILE_ORBITALS_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <votca/tools/matrix.h>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

#include "basis_set.h"
#include "global.h"

namespace votca { namespace moo {

using namespace std;
using std::runtime_error;
using namespace boost;

class orb
{
private :
    
    vector<string> bs;
    double **psi;
    unsigned int NBasis; // the number of basis sets
    int read_orb_gauss(const char *);
    int read_orb_gamess(const char *);
    vector < pair<int, int> > _basis_on_atom;
    basis_set * _basis;
    vector<double> evl;

    double parsenrg(string & line);
public:
    
    static int (orb::*read_orb)(const char *) ;
    orb (): psi(NULL), NBasis(0) {
	//psi=0; Basis=0;
    }

    void assign_basis_set ( basis_set * a, vector< pair <int, int> > b){
        _basis = a;
        vector < pair <int, int> >::iterator it;
        #ifdef DEBUG
        cout << "assigning pairs for so many atoms:" << b.size() <<endl;
        #endif
        for ( it = b.begin(); it != b.end() ; ++it){ 
            pair <int , int > beg_end(it->first, it->second);
            _basis_on_atom.push_back(beg_end);
        }
    }
    
    ~orb(){
	clear();
    }


    void clear(){
    //    cout << "call destructor for orb" <<endl;
    	if (NBasis != 0 ){
	    delete [] psi[0];
	    delete [] psi;
	}
	NBasis =0;
	psi = 0;
	bs.clear();
        evl.clear();
        _basis_on_atom.clear();
    }


    inline const unsigned int & getNBasis() const {
    	return NBasis;
    }

    inline double* getorb(const int & i) const {
    	return psi[i];
    }

    inline double* getevl()  {
    	return &evl[0];
    }

    void cp_orb( orb const  &A){
	for (unsigned int i=0; i < NBasis;i++) {
		for (unsigned int j =0; j < NBasis ; j++){
			psi[i][j] = A.psi[i][j];
		}
	}	
    }

    void cp_orb( orb const &A, int const & i){
	for (unsigned int j =0; j < NBasis ; j++){
	    psi[i][j] = A.psi[i][j];
	}
    }

    void init_orbitals (string * basis, const int & N, const char * namefile );

    void reorder_for_libint();
    
    void set_read_orb_gam();

    void set_read_orb_gauss();
    
    
    void init_orbitals (const orb & orb1 ){
	NBasis = orb1.NBasis;
	bs.resize(NBasis);
	psi = new double* [NBasis];
	psi[0] = new double [NBasis * NBasis];
	bs[0]=orb1.bs[0];
	for ( unsigned int i = 1 ; i < NBasis ; i ++){
		bs[i] = orb1.bs[i];
		psi[i] = psi[i-1] + NBasis;
	}
        for ( unsigned int i =0 ; i< NBasis*NBasis ; ++i) {
            *(psi[0]+i) = *(orb1.psi[0]+i);
        }
        evl.resize(NBasis);
    }
    
    
   
    /*// this function will take the orbitals from two orbitals and generate a guess for the dimer
    *	 all orbitals from both molecules upto homo-1 will be listed first.  
    */
    void dimerise_orbs(const orb & A, const orb & B, const int &elA, const int &elB) ;
    ///this function will print out an UHF guess wavefunction suitable for gaussian usage. It will print
    ///the nel_A/2-1 orbitals (doubly occupied, then the Nbasis -> NBasis+nel_B/2-1 orbitals, then the 
    ///nel_A orbital, then the NBAsis_A+nel_B orbital and finally the rest
    void print_uhf_g03(const int & nel_A, const int & nel_B, const unsigned int & NBasis_A , const unsigned int NBasis_B);

    void print_g03(string & name, string mode= string ("w"));



    void rot_orbs(const vector <unsigned int>& orbs, unsigned int* i, double* psi2, const matrix& rot);
    void rot_orb ( const unsigned int &, const double [3][3]);
    void rot_orb( const double [3][3]);

 
    void rot_orb ( const unsigned int &, const matrix &);
    void rot_orb ( const unsigned int &, unsigned int *, const matrix &);
    void rot_orb ( const unsigned int &, unsigned int *, double *, const matrix &);
    void rot_orb( const matrix &);
    
    void rotate_someatoms(vector<int> , matrix  *, 
            double * , const int &) ;
            
    //void rotate_someatoms(vector<int>, matrix*, double*, const vector <unsigned int>&);
    
    /// removes all orbitals except the a_th
    /*void strip_orbitals( const int & a){
        double * psinew = new double [ NBasis];
        for (int j =0; j < NBasis ; j++){
	    psinew[j] = psi[a][j];
	}
        delete [] psi[0];
        delete [] psi;
        psi = new double *[1];
        psi[0] = psinew;
    }*/
    
    /// removes all orbitals except the relevant transorbs
    void strip_orbitals (const vector < int>& a);
    
    void init_orbitals_stripped(const orb & orb1, const int& nrorbs );
};

}}

#endif 
