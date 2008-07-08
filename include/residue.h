// 
// File:   residue.h
// Author: ruehle
//
// Created on April 10, 2007, 12:01 PM
//

#ifndef _residue_H
#define	_residue_H

#include <string>

using namespace std;
    
/**
    \brief class for a residue
 
    The Residue class describes a residue. When reading in the atoms, all beads belong to a residue. Later on, the molecules
    can be organized into molecules based on their residue.

*/
class Residue
{
public:
    /// constructor
    Residue(int id, const string &name)
        : _id(id), _name(name)
    {}
    
    /// get the name of the residue
    const string &getName();

    /// get the name of the residue
    const int &getId() const { return _id; }

    private:
    int _id;
    string _name;
    
};

inline const string &Residue::getName()
{
    return _name;
}

#endif	/* _residue_H */

