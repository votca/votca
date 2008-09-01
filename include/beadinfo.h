// 
// File:   beadinfo.h
// Author: ruehle
//
// Created on April 5, 2007, 11:17 AM
//


#ifndef _bead_H
#define	_bead_H

#include <string>
#include <tools/types.h>

using namespace std;

/**
    \brief information about a bead
 
    The CBeadInfo class describes an atom or a coarse grained bead. It stores information like the id, the name, the mass, the
    charge and the residue it belongs to. The coordinates are stored in the configuration class.

*/
class BeadInfo
{
public:
    /// constructur
    BeadInfo(int id, byte_t symmetry, string name, int resnr, double m, double q)
        : _id(id), _symmetry(symmetry), _name(name), _resnr(resnr), _m(m), _q(q)
    { }
    
    /// get the id of the bead
    const int &getId() const { return _id; }
    /// get the name of the bead
    const string &getName() const { return _name; }
    
    /// get the residu number of the bead
    const int &getResnr() const { return _resnr; }
    /// get the mass of the bead
    const double &getM() const { return _m; }
    /// get the charge of the bead
    const double &getQ() const { return _q; }
    
    /// get the symmetry of the bead
    /// 1: sphere
    /// 2: ellipsoid with 2 equal axis
    /// 3: ellispoid with 3 different axis
    const byte_t getSymmetry() const { return _symmetry; }

private:
    int _id;
    byte_t _symmetry;
    string _name;
    
    int _resnr;
    
    double _m;
    double _q;
};

#endif	/* _beadinfo_H */

