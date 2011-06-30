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
    reorganization energie, internal coordinates and orbitals.

    \todo mod orbitals to be "stripped down" for memory?

*/

//#define DEBUG

#include "mol_and_orb.h"
#include "orbitals.h"
#include <string>
#include <stdexcept>

class CrgUnit;

class CrgUnitType{
public:
    
    ///it will be initialised by a list of ID names, molecule specifications and orbital
    ///specifications
    
    ~CrgUnitType();
   
    /// reorganization energy of the crg unit type
    const double & getReorg () const { return _reorg; }
    /// radius of chrgunit as sphere for lambda out

    /// site energy of the crg unit type
    const double & getEnergy () const { return _energy; }
    
    const vector <unsigned int>& GetTransOrbs() const;
    
    const unsigned int & getId() const { return _id; }
    
    string & GetMolName();   
    
    const string & GetName() const { return _name; }
    
    // TODO: mol_and_orb should not be in the interface!
    mol_and_orb & GetCrgUnit();
    // TODO: orb should not be in the interface!
    const orb & GetOrb() const;
    // 
    basis_set & GetBS(){
	return _bs;
    }

    Property *getOptions() { return _options; }
    void setOptions(Property *options) { _options = options; }

private:
    basis_set         _bs;
    mol_and_orb       _intcoords;
    orb               _orbitals; 
    double            _reorg;
    double            _energy;

    Property *_options;
    
    vector <unsigned int> _transorbs;
    unsigned int      _id;
    /// the name of this transport unit
    string            _name;  
    // a list of atoms in the same bead of which atoms correspond to which bead
    vector < vector <int>  > _list_atoms_monomer;
    // a list of centre of mass of the monomer
    vector < vec  >   _list_coms_monomer;
    // a list of rotation matrices to put the internal coordinates of each monomer onto the reference state
    vector < matrix > _list_ors_monomer;

    /// this willA take each bead and move it to positions[i] rotating by the
    /// orientation corresponding to norm and rotate we assume that the pointer
    /// is to a suitable molecule...
    void rotate_each_bead(
            vector < vec >::iterator it_pos , vector < vec >::iterator it_norm,
            vector <vec >::iterator it_plan, mol_and_orb * rotated_molecule );

    // can only be created by JCalc
    CrgUnitType(): _reorg(0.), _energy(0.), _id(-1) {};
    
    // can only be created by JCalc
    CrgUnitType(
            const char * namecoord, const char * nameorb, const char * nameneutr,
            const char * namecrg, string & basisset, const double & reorg,
            const double & energy, const vector < int>& transorbs,
            const unsigned int &id,  string name,
            vector < vector < int > > list_atoms_monomer,
            vector < vector < double > > list_weights_monomer);

    friend class CrgUnit;
    friend class JCalc;
};

inline const vector <unsigned int>& CrgUnitType::GetTransOrbs() const{
    return _transorbs;
}



inline mol_and_orb & CrgUnitType::GetCrgUnit() {
    return _intcoords;
}

inline const orb & CrgUnitType::GetOrb() const {
    return _orbitals;
}

#endif	/* _CRGUNITTYPE_H */

