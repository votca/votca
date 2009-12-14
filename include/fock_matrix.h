#ifndef FILE_FOCK_H
#define FILE_FOCK_H


#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "orbitals.h"
#include "global.h"
#include "qm_molecule.h"



/// this class will contain the fock matrix and pointers to the molecules involved in calculating it
class fock {
    private:

	pair <const mol_and_orb *, const mol_and_orb *> molecules; // an ordered pair of pointers to the molecules
	double *Beta, *Mu; // the choice of paramters

	double **F; 

	void calc_F_el_with_HASH( const double & , double & , const double & , const double & , const double & , const double & , 
		const int & , const int & , const int & , const int &  ) const;
	void calc_F_lean( ) const ;

        static bool _init_HASH;

    public:

	/// the fock matrix to calculate F[i][j] will contain the ith elemnt of the basis on the first molecule and the jth element on the second  
	
	fock () {
		F = NULL;
		molecules.first = NULL;
		molecules.second = NULL;
		
	}

	fock (  const mol_and_orb & A, const mol_and_orb & B ){
		molecules.first  = &A;
		molecules.second = &B;
		int n1, n2;

		n1 = molecules.first  -> getNBasis();
		n2 = molecules.second -> getNBasis();

		F = new double* [n1];
		F[0] = new double [n1*n2];
		for (int i = 1 ; i < n1 ; ++i){
		    F[i] = F[i-1] +  n2;
		}
		
		if ( _init_HASH == false ) {
		    fock::init_HASH();
                    _init_HASH = true;
		}
	}

        void init(  const mol_and_orb & A, const mol_and_orb & B ){
		molecules.first  = &A;
		molecules.second = &B;
		int n1, n2;

		n1 = molecules.first  -> getNBasis();
		n2 = molecules.second -> getNBasis();

		F = new double* [n1];
		F[0] = new double [n1*n2];
		for (int i = 1 ; i < n1 ; ++i){
		    F[i] = F[i-1] +  n2;
		}
		
		if ( _init_HASH == false ) {
		    fock::init_HASH();
                    _init_HASH = true;
		}
                set_zindo_s();
	}
        
	~fock() {
            
        #ifdef DEBUG
        cout << "callgin the fockk classdestructor" << endl;
        #endif
	    clear();
	}

	void clear () {
	    if ( F != NULL ) {
		delete [] F[0];
		delete [] F;
		F = NULL;
	    }
	}

	void set_zindo_1();
	void set_zindo_s();

	vector <double> calcJ( vector<pair <unsigned int, unsigned int> >) const;
	double calcJ( pair <int, int>) const;
	
	int init_HASH();

};

#endif
