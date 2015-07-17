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

#ifndef FILE_MOLECULES_EL
#define FILE_MOLECULES_EL

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip> 
#include <vector>
#include <stdlib.h>
#include <cmath>
#include "orbitals.h"
#include "basis_set.h"
#include "charges.h"
#include "global.h"
#include "units.h"
#include <votca/tools/vec.h>
#include <votca/tools/globals.h>

////#define DEBUG
namespace votca { namespace moo {

using namespace std;

class mol_and_orb
{
private:

    struct atom_type {
	string _type;
	int _lbl;
    };   

    vector <atom_type> atom_labels ; // label of atoms
    vector <vec>       atom_pos;
    unsigned int N; // the number of atoms in the molecule

    vec centre_input;

    orb *orbitals;
    multipoles *_crged;
    multipoles *_neutr;
    vec centre;
    basis_set *_bs;

public:

    int n_el;

    mol_and_orb(): N(0), orbitals(0), _crged(0),
                   _neutr(0), centre(0.,0.,0.), n_el(0) { }

   ~mol_and_orb() {
	atom_pos.clear();
	atom_labels.clear();
    }


    vec getCentre() { return centre; }
    vec getCentreInput() { return centre_input; }
    void setCentre(const vec &pos) { centre = pos; }
    vector< vec > getAtomPos() { return atom_pos; }
    void setAtomPos(vector< vec > &positions) { atom_pos = positions; }
    vector< atom_type > getAtomLabels() { return atom_labels; }
    void setAtomLabels(vector<string> types, vector<int> labels);


    multipoles *getChargesNeutr() { return _neutr; }
    multipoles *getChargesCrged(){ return _crged; }

    void cp_atompos(mol_and_orb & A) {
       atom_pos.clear();	 
    
       vector < vec >:: iterator iter_pos;
       for (iter_pos=(A.atom_pos).begin();
            iter_pos<A.atom_pos.end();
            ++iter_pos){
            atom_pos.push_back(*iter_pos);
       }	   
       N = A.N;
    }

    void cp_atoms(mol_and_orb  & A) {
	atom_labels.clear();	 

	vector < atom_type >:: iterator iter_pos;
        for (iter_pos=A.atom_labels.begin();
             iter_pos<A.atom_labels.end();
             ++iter_pos){
             atom_labels.push_back(*iter_pos);
        }	   
    }

    /*void cp_orb( mol_and_orb const & A, vector <int> a ){
    	for (unsigned  int i=0 ;i< a.size(); i++){
		orbitals -> cp_orb( *(A.orbitals), a[i]) ; 
	}
    }*/
    
    void cp_orb( mol_and_orb const & A, vector <unsigned int> a ) {
    	for (unsigned  int i=0 ;i< a.size(); i++) {
             orbitals -> cp_orb( *(A.orbitals), a[i]) ;
	}
    }

    void cp_orb( mol_and_orb const & A, const int & a ) {
        orbitals -> cp_orb( *(A.orbitals), a);
    }

    void cp_crg ( const mol_and_orb & A ) {
        if (A._crged == NULL) return;
        if (A._neutr == NULL) return;
        if (_crged == 0 ) _crged = new multipoles;
        if (_neutr == 0 ) _neutr = new multipoles;
        
    	_crged -> cp_mpls ( * (A._crged));
        _neutr -> cp_mpls ( * (A._neutr));
    }

    void rotate(const matrix &M ) {
	vector < vec >:: iterator iter_pos;
	for (iter_pos = atom_pos.begin();
             iter_pos <  atom_pos.end();
             ++iter_pos){
             *iter_pos = M * (*iter_pos);
	}
	centre = M * centre;
    }
    
    void rot_orb(const double rot[3][3], vector <unsigned int> a ) {
    	for (unsigned int i=0; i< a.size();i++) {
             orbitals -> rot_orb(a[i], rot);
	}
    }

    void rot_orb( const double rot[3][3], const unsigned int & i) {
	orbitals -> rot_orb(i, rot);	    
    }

    
    void rot_orb( const matrix &rot, vector <unsigned int> a ) {
    	for (unsigned int i=0; i< a.size();i++) {
	     orbitals -> rot_orb(a[i], rot);
	}
    }

    void rot_orb( const matrix &rot, const unsigned int & i) {
	orbitals -> rot_orb(i, rot);	    
    }
    
    void shift(double a[3]) {
        vec A(a);
        shift(A);
    }

    void shift(const vec &a ) {
	vector < vec >:: iterator iter_pos;
	for (iter_pos = atom_pos.begin();
             iter_pos <  atom_pos.end();
             ++iter_pos) {
            *iter_pos += a;
	}
    }    
    
    void print(ostream & out) {
	out.setf(ios::fixed);
     	for (unsigned int i=0; i < N ;i++ ) {
	    out << atom_labels[i]._type << '\t' 
                << unit<bohr,A>::to(atom_pos[i].getX()) << '\t'
                << unit<bohr,A>::to(atom_pos[i].getY()) << '\t'
                << unit<bohr,A>::to(atom_pos[i].getZ()) << '\t'
                << '\n' ;
	}
    }
     
    void print() {
        print(cout);
    }
		
    inline int getlbl (const int & i ) const {
	return (atom_labels[i])._lbl;
    }
		
    inline string gettype (const int & i ) const {
	return (atom_labels[i])._type;
    }

    inline const unsigned int & getN() const {
    	return N;
    }

    inline const unsigned int & getNBasis() const {
    	return orbitals->getNBasis();
    }

    orb * getorbs() const {
    	return orbitals;
    }
    
    /*inline double * getorb(const int & i) const {
    	return orbitals->getorb(i);
    }*/
 
    inline double * getorb(const unsigned int & i) const {
    	return orbitals->getorb(i);
    }

    inline const double & getx(const int & i) const {
    	return (atom_pos[i]).getX();
    }

    inline const double & gety(const int & i) const {
    	return (atom_pos[i]).getY();
    }

     inline const double & getz( const int & i) const {
    	return (atom_pos[i]).getZ();
    }
     
     inline const vec & GetPos ( const int & i ) const {
	 return atom_pos[i];
     }
     
     
     inline const int  & get_nel_at(const int & i) const {
         return _bs->get_nel_at(atom_labels[i]._lbl);
     }
     
     inline const int  & get_nbasis_at(const int & i) const {
         return _bs->get_nbasis_at(atom_labels[i]._lbl);
     }

     basis_set * get_basis() const{
	 return _bs;
     }

    /*const double & get_mpls(const int & i ) const{
     	return _charge->get_mpl(i);
     }*/
     
     ///makes a list of start orbital on an atom and number if basis sets on that atom
     void assign_orb ( orb * orb1){
         orbitals = orb1;
         //determine the vector of pairs.
         vector < pair < int, int> > basis_on_atom;
         vector <atom_type>::iterator it_at_type;
         int beg = 0;
         //int nbasis;
         n_el = 0;
         for (it_at_type = atom_labels.begin();
              it_at_type != atom_labels.end();
              ++it_at_type){

             pair <int, int> basis_fo; 
             basis_fo.first = beg;
             basis_fo.second = _bs -> get_nbasis_at(it_at_type -> _lbl);
             basis_on_atom.push_back(basis_fo);
             beg += basis_fo.second;

             n_el +=_bs->get_nel_at(it_at_type -> _lbl);
         }
         orb1 -> assign_basis_set(_bs, basis_on_atom);
     }
     
     void assign_charge (multipoles * crg1, multipoles * crg2) {
	 _crged = crg1;
         _neutr = crg2;
     }

   int init(const char * );

   int init( char * a) {
	const char * A = a;
	init (A);
	return 0;
    }

   const int &getNEl() { return n_el; }
   

   int init_orbitals (orb &,  const char * );
   int init_charges ( const char *, const char *);
   
   void define_bs ( basis_set & bs) {
       _bs = &bs;
   }
   void define_bs ( basis_set *bs) {
       _bs = bs;
   }
   void set_basis_set(const string & a){
       _bs->set_basis_set(a);
   }

    double V_mpls_neutr(vec &a){ //returns V for potential
	double V=0.0; 
	vec dist;
	for (unsigned int i=0;i<N;i++) {
	     dist = a-atom_pos[i];
	     double r = abs(dist);
	     V +=  _neutr->get_mpl(i) / r;
	}
	return V;
    }
    
    double V_mpls_crg(vec &a){ //returns V for potential
	double V=0.0; 
	vec dist;
	for (unsigned int i=0;i<N;i++) {
	     dist = a-atom_pos[i];
	     double r = abs(dist);
	     V +=  _crged->get_mpl(i) /  ( r );
	}
	return V;
    }
    
    double V_nrg_neutr_neutr(mol_and_orb & a){
        double V=0;
        for (unsigned int i=0; i<N;i++){
            V += _neutr->get_mpl(i) * a.V_mpls_neutr(atom_pos[i]);
        }
        return V;
    }
    
    double V_nrg_crg_neutr(mol_and_orb & a){
        double V=0;
        for (unsigned int i=0; i<N;i++){
            V += _crged->get_mpl(i) * a.V_mpls_neutr(atom_pos[i]);
        }
        return V;
    }
    
    void write_pdb(string,string, const int &);

    void rotate_someatoms(vector <int> list_at,
                          const matrix  & M,
                          const vec & a,
                          const vec & com,
                          mol_and_orb * mol2) {
        vector <int>::iterator it_at ;
        
        if (globals::verbose) {
                cout << "translating " << list_at.size() 
                     << " atoms by distance: " << a << endl <<
                        " rotating by " << endl 
                     << M.get(0,0) << " " << M.get(0,1) << " " << M.get(0,2)  << endl
                     << M.get(1,0) << " " << M.get(1,1) << " " << M.get(1,2)  << endl
                     << M.get(2,0) << " " << M.get(2,1) << " " << M.get(2,2)  << endl;   
        }
        
            #ifdef DEBUG
                 cout << "about to rotate so many atoms: "
                 << list_at.size() << " by this distance: "
                 << a << endl;
            #endif

        for (it_at = list_at.begin(); it_at != list_at.end(); ++it_at) {

            #ifdef DEBUG
                 cout << "rotating atom:" << *it_at <<endl;
            #endif           

            mol2->atom_pos[*it_at].setX( a.getX() + M.get(0,0) * (atom_pos[*it_at].getX() - com.getX() ) + M.get(0,1) * (atom_pos[*it_at].getY() - com.getY() ) + M.get(0,2) * (atom_pos[*it_at].getZ() - com.getZ() ));
            mol2->atom_pos[*it_at].setY( a.getY() + M.get(1,0) * (atom_pos[*it_at].getX() - com.getX() ) + M.get(1,1) * (atom_pos[*it_at].getY() - com.getY() ) + M.get(1,2) * (atom_pos[*it_at].getZ() - com.getZ() ));
            mol2->atom_pos[*it_at].setZ( a.getZ() + M.get(2,0) * (atom_pos[*it_at].getX() - com.getX() ) + M.get(2,1) * (atom_pos[*it_at].getY() - com.getY() ) + M.get(2,2) * (atom_pos[*it_at].getZ() - com.getZ() ));
        }
    }


    void rotate_someatoms_ctp(vector <int> list_at,
                              const matrix  & M,
                              const vec & a,
                              const vec & com,
                              mol_and_orb * mol2) {
        vector <int>::iterator it_at ;

        for (it_at = list_at.begin(); it_at != list_at.end(); ++it_at) {
            mol2->atom_pos[*it_at] = M * ( mol2->atom_pos[*it_at] - com ) + a;
        }
    }
};

}}

#endif //FILE_MOLECULES_EL
