/// \addtogroup csg
///@{
// 
// File:   configuration.cc
// Author: ruehle
//
// Created on April 5, 2007, 12:30 PM
//

#include "topology.h"
#include "configuration.h"

Configuration::Configuration(Topology *top)
        : _topology(top)
{
    Cleanup();
}

Configuration::~Configuration()
{
    Cleanup();
}

void Configuration::Initialize()
{
    _nBeads = _topology->BeadCount();
}

void Configuration::Cleanup()
{
    _bead_positions.clear();
    _bead_velocities.clear();
    _bead_forces.clear();
    _bead_u.clear();
    _bead_v.clear();
    _bead_w.clear();

    _bPos=_bVel=_bF=_bU=_bV=_bW=false;
}

/// \todo check weather it really is pbc
vec Configuration::getDist(int bead1, int bead2) const
{
    vec r_tp, r_dp, r_sp, r_ij;
    vec r_i = _bead_positions[bead1];
    vec r_j = _bead_positions[bead2];
    vec a = _box.getCol(0); vec b = _box.getCol(1); vec c = _box.getCol(2);
    r_tp = r_j - r_i;
    r_dp = r_tp - c*round(r_tp.getZ()/c.getZ());  
    r_sp = r_dp - b*round(r_dp.getY()/b.getY());
    r_ij = r_sp - a*round(r_sp.getX()/a.getX());
    return r_ij;
}

double Configuration::BoxVolume()
{
    vec a = _box.getCol(0); vec b = _box.getCol(1); vec c = _box.getCol(2);
    return (a^b)*c;
}
/// @}
