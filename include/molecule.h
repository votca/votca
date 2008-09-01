// 
// File:   molecule.h
// Author: ruehle
//
// Created on April 17, 2007, 2:14 PM
//

#ifndef _molecule_H
#define	_molecule_H

#include "configuration.h"
#include "moleculeinfo.h"

/**
    \brief molecule - connection between moleculeinfo and a configuration

    The CMolecule class connects a molecule info with a configuration to represent the whole molecule.
*/
class Molecule
{
public:
    /// constructor
    Molecule(Configuration &conf, MoleculeInfo &info)
        : _conf(&conf), _info(&info) {}
        
    /// get the configuration which corresponds to the molecule
    Configuration *getConfiguration() { return _conf; }
    /// get the molecule info which stores the information about the molecule
    MoleculeInfo *getMoleculeInfo() { return _info; }
    
    /// wrapper to easily read the position of a bead in the molecule
    const vec &getBeadPos(const int &bead) const;
    /// wrapper to get direct access (read/write) to the position of a bead in the molecule
    vec &BeadPos(const int &bead);
    
    const vec &getBeadF(const int &bead) const;
    vec &BeadF(const int &bead);
    
    const vec &getBeadVel(const int &bead) const;
    vec &BeadVel(const int &bead);
    
    vec &BeadU(const int &bead);
    vec &BeadV(const int &bead);
    vec &BeadW(const int &bead);
    
private:
    Configuration *_conf;
    MoleculeInfo *_info;
};

inline const vec &Molecule::getBeadPos(const int &bead) const
{
    return _conf->Pos(_info->getBeadId(bead));
}

inline vec &Molecule::BeadPos(const int &bead)
{
    return _conf->Pos(_info->getBeadId(bead));
}

inline const vec &Molecule::getBeadF(const int &bead) const
{
    return _conf->F(_info->getBeadId(bead));
}

inline const vec &Molecule::getBeadVel(const int &bead) const
{
    return _conf->Vel(_info->getBeadId(bead));
}

inline vec &Molecule::BeadVel(const int &bead)
{
    return _conf->Vel(_info->getBeadId(bead));
}

inline vec &Molecule::BeadF(const int &bead)
{
    return _conf->F(_info->getBeadId(bead));
}


inline vec &Molecule::BeadU(const int &bead)
{
    return _conf->U(_info->getBeadId(bead));
}

inline vec &Molecule::BeadV(const int &bead)
{
    return _conf->V(_info->getBeadId(bead));
}

inline vec &Molecule::BeadW(const int &bead)
{
    return _conf->W(_info->getBeadId(bead));
}

#endif	/* _molecule_H */

