/// \addtogroup csg
///@{
// 
// File:   moleculeinfo.h
// Author: ruehle
//
// Created on April 12, 2007, 4:07 PM
//

#ifndef _moleculeinfo_H
#define	_moleculeinfo_H

#include <vector>
#include <map>
#include <string>
#include <assert.h>

using namespace std;
    
/**
    \brief information about molecules

    The CMoleculeInfo class stores which beads belong to a molecule.
    The organization of beads into molecules is needed for the CG mapping.

    \todo sort atoms in molecule

*/
class MoleculeInfo
{
public:        
    /// constructor
    MoleculeInfo(int id, string name)
        : _id(id), _name(name)
    {}
    
    /// get the molecule ID
    int getId() const { return _id; }
    
    /// get the name of the molecule
    const string &getName() const { return _name; }
    
    /// set the name of the molecule
    void setName(const string &name) {  _name=name; }
    
    /// Add a bead to the molecule
    void AddBead(int bead, const string &name);
    /// get the id of a bead in the molecule
    int getBeadId(int bead) { return _beads[bead]; }
    int getBeadIdByName(const string &name) { return _beads[getBeadByName(name)]; }
    
    /// get the number of beads in the molecule
    int BeadCount() const { return _beads.size(); }
    
    /// find a bead by it's name
    int getBeadByName(const string &name);
    string getBeadName(int bead) {return _bead_names[bead]; }
    
private:
    // maps a name to a bead id
    map<string, int> _beadmap;
    
    // id of the molecules
    int _id;
    
    // name of the molecule
    string _name;
    // the beads in the molecule
    vector<int> _beads;
    vector<string> _bead_names;
};

#endif	/* _moleculeinfo_H */

/// @}
