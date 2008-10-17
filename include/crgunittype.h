// 
// File:   crgunittype.h
// Author: james
//
// Created on 27 June 2007, 09:21
//

#ifndef _CRGUNITTYPE_H
#define	_CRGUNITTYPE_H
/**
    \brief information about a CrgUnitTpe
 
    The CrgUnitType class describes a known charge unit type. It stores the ID,
 * reorganization energie, internal coordinates and orbitals. 

    \todo mod orbitals to be "stripped down" for memory?

*/

//#define DEBUG

#include "qm_molecule.h"
#include "orbitals.h"
#include <libxml/parser.h>
#include <string>

inline void get_orient( const vec & a, const vec & b, const vec & c, matrix & cg ){
    
    cg.set(0,0,  a.getX());
    cg.set(0,1,  b.getX());
    cg.set(0,2,  c.getX());
    cg.set(1,0,  a.getY());
    cg.set(1,1,  b.getY());
    cg.set(1,2,  c.getY());
    cg.set(2,0,  a.getZ());
    cg.set(2,1,  b.getZ());
    cg.set(2,2,  c.getZ());
}

class CrgUnitType{
public:
    
    ///it will be initialised by a list of ID names, molecule specifications and orbital
    ///specifications
    
    CrgUnitType(){};
    
    ~CrgUnitType(){
        #ifdef DEBUG
        cout << "callgin the crgunitype destructor" << endl;
        #endif
        vector < vector <int>  >::iterator  it_list_atoms_monomer = _list_atoms_monomer.begin();
        for ( ; it_list_atoms_monomer !=  _list_atoms_monomer.end() ; ++ it_list_atoms_monomer){
            it_list_atoms_monomer->clear();
        }
        _list_atoms_monomer.clear();
        _list_coms_monomer.clear();
        _list_ors_monomer.clear();
    }

    CrgUnitType (char * namecoord, char * nameorb, char * nameneutr, char * namecrg, const double & reorg, 
            const double & energy, const unsigned int &transorb, const unsigned int &id, 
            string molname, string name, vector < vector < int > > list_atoms_monomer ){
        _intcoords.define_bs(_indo);
        _intcoords.init(namecoord);
        _intcoords.init_orbitals( _orbitals, nameorb);
        _intcoords.init_charges(nameneutr, namecrg);
        #ifdef DEBUG
        cout << "sample of the trans: " << transorb << " orbitals: " << (_orbitals.getorb(transorb))[4] <<endl;
        #endif
        _orbitals.strip_orbitals(transorb);
        #ifdef DEBUG
        cout << "sample of the stripped orbitals: " << (_orbitals.getorb(0))[4] <<endl;
        #endif
        _reorg    = reorg;
        _energy   = energy;
        _transorb = 0;
        _id       = id;
        _molname  = molname;
        _name     = name;
        vector < vector < int > >  :: iterator it_mon;
        int count=0;
        for (it_mon= list_atoms_monomer.begin() ; it_mon < list_atoms_monomer.end() ; ++it_mon ){
            vector <int> list_mon;
            _list_atoms_monomer.push_back(list_mon);
            vector <int> ::iterator  it_at;
            vec com(0.,0.,0.);
            matrix m(0.);
            matrix orient;
            vec xprime, yprime, zprime;
            for(it_at= it_mon->begin(); it_at != it_mon->end(); ++it_at){
                com = com + _intcoords.GetPos(*it_at);
                (_list_atoms_monomer[count]).push_back(*it_at);
            }

            com /= it_mon->size();
            _list_coms_monomer.push_back(com);
            
            if (it_mon-> size() > 3 ){
                for (it_at = it_mon->begin(); it_at != it_mon->end(); ++it_at){
                    vec v = _intcoords.GetPos(*it_at);
                    #ifdef DEBUG
                    cout << "Adding atom " << *it_at << " to monomer " <<count <<endl;
                    #endif
                    m[0][0] += v.getX()*v.getX();
                    m[0][1] += v.getX()*v.getY();
                    m[0][2] += v.getX()*v.getZ();        
                    m[1][1] += v.getY()*v.getY();
                    m[1][2] += v.getY()*v.getZ();
                    m[2][2] += v.getZ()*v.getZ();        
                }

                m[1][0] = m[0][1];
                m[2][0] = m[0][2];
                m[2][1] = m[1][2];

                matrix::eigensystem_t es;
                m.SolveEigensystem(es);
                zprime = es.eigenvecs[0];


                xprime = _intcoords.GetPos((list_atoms_monomer[count])[1]) - _intcoords.GetPos((list_atoms_monomer[count])[0] );
                xprime.normalize(); //these 2 lines useless
                zprime.normalize();
                #ifdef DEBUG 
                cout << "xaxis : " << xprime <<endl;
                #endif
                vec w = _intcoords.GetPos((list_atoms_monomer[count])[2]) - _intcoords.GetPos((list_atoms_monomer[count])[0] );
                w.normalize();
                #ifdef DEBUG
                cout << "CN bond: " <<  w <<endl;
                #endif
                if((xprime^w)*zprime < 0) {
                    zprime=vec(0.,0.,0.)-zprime;
                    #ifdef DEBUB
                    cout << "flipping the zaxis" <<endl;
                    #endif
                }

                vec yprime = (zprime) ^ (xprime);
                xprime = yprime ^ (zprime);
                xprime.normalize();
                yprime.normalize();
                zprime.normalize();

                
                get_orient(xprime, yprime, zprime, orient);
                orient.Transpose();
            }
            else {
                
                orient[0][0] = 1.0;
                orient[0][1] = 0.;
                orient[0][2] = 0.;
                orient[1][0] = 0.;
                orient[1][1] = 1.0;
                orient[1][2] = 0.;
                orient[2][0] = 0.;
                orient[2][1] = 0.;
                orient[2][2] = 1.0;
            }
            _list_ors_monomer.push_back(orient);
            
            #ifdef DEBUG 
            cout << "Making charge unit of type: "  << _id << " initing monomer " << count << " with normal " 
                    << zprime << " and plane: " << xprime << endl;
            #endif
            count++;
        }
        
        _intcoords.assign_orb(&_orbitals);
    }
    
    inline const double & GetReorg () const{
        return _reorg;
    }
    
    inline const double & GetEnergy () const {
        return _energy;
    }
    
    inline const unsigned int & GetTransOrb() const{
        return _transorb;
    }
            
    inline const unsigned int & GetId() const{
        return _id;
    }
    
    string & GetMolName(){
        return _molname;
    }
    
    
    string & GetName(){
        return _name;
    }
    
    inline mol_and_orb & GetCrgUnit() {
        return _intcoords;
    }
    
    inline const orb & GetOrb() const{
        return _orbitals;
    }
    /// this willA take each bead and move it to positions[i] rotating by the orientation corresponding to norm and rotate
    /// we assume that the pointer is to a suitable molecule... 
    void rotate_each_bead( vector < vec >::iterator it_pos , vector < vec >::iterator it_norm, 
            vector <vec >::iterator it_plan, mol_and_orb * rotated_molecule ){
        vector < vector <int> >::iterator it_monomers; 
        vector < vec >::iterator it_coms = _list_coms_monomer.begin();
        vector <matrix>::iterator it_ors = _list_ors_monomer.begin();
                
        int count =0;
        vec old_norm;
        for (it_monomers = _list_atoms_monomer.begin(); it_monomers != _list_atoms_monomer.end(); 
        ++it_monomers , ++it_pos, ++it_norm, ++it_plan, ++it_coms, ++it_ors){
            
            
            // we need to define the rotation matrix for the beads
            vec zprime = (*it_norm);
            vec yprime = zprime ^ (*it_plan);
            vec xprime = yprime ^ (zprime);
            
            xprime = xprime.normalize();
            #ifdef DEBUG 
            if (abs(yprime) < 1E-6 ){
                cout << "Errorm the inplane vector is parallel to the nomal one!" <<endl;
            }
            #endif
            yprime = yprime.normalize();
           
            vec zprime_orb = zprime;
            if (count > 0){
                // check if the orbitals need "flipping"
                if ( old_norm * ( zprime) > 0. ) zprime_orb=vec(0.,0.,0.)-zprime;
            }
            
            matrix orient_pos;
            matrix orient_orb;
            
            get_orient(xprime, yprime, zprime, orient_pos);
            get_orient(xprime, yprime, zprime_orb, orient_orb);
            
            
            matrix Orient_Pos = orient_pos * (*it_ors);
            matrix Orient_Orb = orient_orb * (*it_ors);
            
            #ifdef DEBUG
            cout << "NUMBER OF atoms in the the monemer: " << it_monomers->size() << "displacement" << *it_pos <<  
                    " internal cebtre of nassL " << *it_coms << endl;
            cout << " type of crgunittype: " << _id<<endl;
            cout << "normal vector: " << zprime << " and for orbitals: " << zprime_orb<< " and plane vector: " << xprime <<endl;
            #endif
             _intcoords.rotate_someatoms(*it_monomers, Orient_Pos, *it_pos, *it_coms, rotated_molecule);
            _orbitals.rotate_someatoms(*it_monomers, &Orient_Orb, rotated_molecule->getorb(_transorb), _transorb );
            
            old_norm = zprime;
            count++;
        }
    }
    
        
private:
    basis_set         _indo;
    mol_and_orb       _intcoords;
    orb               _orbitals; 
    double            _reorg;
    double            _energy;
    unsigned int      _transorb;
    unsigned int      _id;
    ///the molecule name that this charge tranpost unit belongs to
    string            _molname;
    /// the name of this transport unit
    string            _name;  
    // a list of atoms in the same bead of which atoms correspond to which bead
    vector < vector <int>  > _list_atoms_monomer;
    // a list of centre of mass of the monomer
    vector < vec  >   _list_coms_monomer;
    // a list of rotation matrices to put the internal coordinates of each monomer onto the reference state
    vector < matrix > _list_ors_monomer;
    
};



#endif	/* _CRGUNITTYPE_H */

