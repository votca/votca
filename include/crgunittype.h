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

#include "mol_and_orb.h"
#include "orbitals.h"
#include <libxml/parser.h>
#include <string>
#include <stdexcept>

class CrgUnitType{
public:
    
    ///it will be initialised by a list of ID names, molecule specifications and orbital
    ///specifications
    
    CrgUnitType(){};
    
    ~CrgUnitType();

    CrgUnitType (
            char * namecoord, char * nameorb, char * nameneutr,
            char * namecrg, const double & reorg, 
            const double & energy, const vector <unsigned int>& transorbs,
            const unsigned int &id, string molname, string name,
            vector < vector < int > > list_atoms_monomer,
            vector < vector < double > > list_weights_monomer );

    const double & GetReorg () const;
    const double & GetEnergy () const;
    
    //const unsigned int & GetTransOrb() const;
      
    const vector <unsigned int>& GetTransOrbs() const;
    
    const unsigned int & GetId() const;
    
    string & GetMolName();   
    
    string & GetName();
    
    mol_and_orb & GetCrgUnit();
    
    const orb & GetOrb() const;
    
    /// this willA take each bead and move it to positions[i] rotating by the
    /// orientation corresponding to norm and rotate we assume that the pointer
    /// is to a suitable molecule... 
    void rotate_each_bead(
            vector < vec >::iterator it_pos , vector < vec >::iterator it_norm,
            vector <vec >::iterator it_plan, mol_and_orb * rotated_molecule );
        
private:
    basis_set         _indo;
    mol_and_orb       _intcoords;
    orb               _orbitals; 
    double            _reorg;
    double            _energy;
    //unsigned int      _transorb;
    vector <unsigned int> _transorbs;
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

inline const double & CrgUnitType::GetReorg() const {
    return _reorg;
}

inline const double & CrgUnitType::GetEnergy() const {
    return _energy;
}

/*inline const unsigned int & CrgUnitType::GetTransOrb() const {
    return _transorb;
}*/

inline const vector <unsigned int>& CrgUnitType::GetTransOrbs() const{
    return _transorbs;
}

inline const unsigned int & CrgUnitType::GetId() const {
    return _id;
}

inline string & CrgUnitType::GetMolName() {
    return _molname;
}

inline string &CrgUnitType::GetName() {
    return _name;
}

inline mol_and_orb & CrgUnitType::GetCrgUnit() {
    return _intcoords;
}

inline const orb & CrgUnitType::GetOrb() const {
    return _orbitals;
}

#endif	/* _CRGUNITTYPE_H */

