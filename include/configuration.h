// 
// File:   configuration.h
// Author: ruehle
//
// Created on April 5, 2007, 11:48 AM
//

#ifndef _configuration_H
#define	_configuration_H

#include <vector>
#include <assert.h>
#include <tools/matrix.h>
#include <tools/vec.h>

using namespace std;

class Topology;

/**
    \brief configuration of a system

    This class stores the configuration of the system like bead positions, orientations, ...

    \todo think about weather bead struct which combines pos, uvw, f into one class and store array of that
 
*/
class Configuration
{
public:
    /// constructor
    Configuration(Topology *top);
    ~Configuration();
    
    /// initialize the configuration based on the topology it was created
    void Initialize();
    /// cleanup
    void Cleanup();
    
    /// get the topology
    Topology *getTopology() { return _topology; }
    
    /// set position of a bead
    void setPos(const int &bead, const vec &r);
    /// get the position of a bead
    const vec &getPos(const int &bead) const;

    /// set velocity of a bead
    void setVel(const int &bead, const vec &r);
    /// get the velocity of a bead
    const vec &getVel(const int &bead) const;
    
    /// set orientation U of a bead
    void setU(const int &bead, const vec &u);
    /// get orientation U of a bead
    const vec &getU(const int &bead) const;

    /// set orientation V of a bead
    void setV(const int &bead, const vec &v);
    /// get orientation V of a bead
    const vec &getV(const int &bead) const;

    /// set orientation W of a bead
    void setW(const int &bead, const vec &w);
    /// get orientation W of a bead
    const vec &getW(const int &bead) const;
    
    
    /// direct access (read/write) to the position of a bead
    vec &Pos(const int &bead) { return _bead_positions[bead]; }
    /// direct access (read/write) to the velocities of a bead
    vec &Vel(const int &bead) { return _bead_velocities[bead]; }
    /// direct access (read/write) to the orientation U of a bead
    vec &U(const int &bead) { return _bead_u[bead]; }
    /// direct access (read/write) to the orientation U of a bead
    vec &V(const int &bead) { return _bead_v[bead]; }
    /// direct access (read/write) to the orientation U of a bead
    vec &W(const int &bead) { return _bead_w[bead]; }
    /// direct access (read/write) to the force of a bead
    vec &F(const int &bead) { return _bead_forces[bead]; }
 
    /// set force acting on a bead
    void setF(const int &bead, const vec &F);
    /// get the force acting on a bead
    const vec &getF(const int &bead) const;
    
    /// set the simulation box
    void setBox(const matrix &box) { _box = box; };
    /// get the simulation box
    const matrix &getBox() { return _box; };

    /// set the time of current frame
    void setTime(real t) { _time = t; };
    /// get the time of current frame
    real getTime() { return _time; };
    
    /// set the current step
    void setStep(int s) { _step = s; };
    /// get the current step
    int getStep() { return _step; };

    /// get the smallest distance between two particles with correct treatment of pbc
    vec getDist(int bead1, int bead2) const;
    
    /// calculates the vox volume
    double BoxVolume();
    
    /// does this configuration store positions?
    bool HasPos() {return _bPos; }
    /// does this configuration store velocities?
    bool HasVel() {return _bVel; }
    /// does this configuration store forces?
    bool HasF() {return _bF; }
    /// does this configuration store u-orientations?
    bool HasU() {return _bU; }
    /// does this configuration store v-orientations?
    bool HasV() {return _bV; }
    /// does this configuration store w-orientations?
    bool HasW() {return _bW; }
    
    void HasPos(bool b);
    void HasVel(bool b);
    void HasF(bool b);
    void HasU(bool b);
    void HasV(bool b);
    void HasW(bool b);
    
private:
    Topology *_topology;
    matrix _box;
    real _time;
    int _step;
    int _nBeads;
    
    vector<vec> _bead_positions;
    vector<vec> _bead_velocities;
    vector<vec> _bead_forces;
    vector<vec> _bead_u;
    vector<vec> _bead_v;
    vector<vec> _bead_w;
    
    bool _bPos;
    bool _bVel;
    bool _bU;
    bool _bV;
    bool _bW;
    bool _bF;
};

inline void Configuration::setPos(const int &bead, const vec &r)
{
    assert(_bPos);
    _bead_positions[bead] = r;
}

inline const vec &Configuration::getPos(const int &bead) const
{
    assert(_bPos);
    return _bead_positions[bead];
}

inline void Configuration::setVel(const int &bead, const vec &r)
{
    assert(_bVel);
    _bead_velocities[bead] = r;
}

inline const vec &Configuration::getVel(const int &bead) const
{
    assert(_bVel);
    return _bead_velocities[bead];
}

inline void Configuration::setU(const int &bead, const vec &u)
{
    assert(_bU);
    _bead_u[bead] = u;
}

inline const vec &Configuration::getU(const int &bead) const
{
    assert(_bU);
    return _bead_u[bead];
}

inline void Configuration::setV(const int &bead, const vec &v)
{
    assert(_bV);
    _bead_v[bead] = v;
}

inline const vec &Configuration::getV(const int &bead) const
{
    assert(_bV);
    return _bead_v[bead];
}

inline void Configuration::setW(const int &bead, const vec &w)
{
    assert(_bW);
    _bead_w[bead] = w;
}

inline const vec &Configuration::getW(const int &bead) const
{
    assert(_bW);
    return _bead_w[bead];
}

inline void Configuration::setF(const int &bead, const vec &F)
{
    assert(_bF);
    _bead_forces[bead] = F;
}

inline const vec &Configuration::getF(const int &bead) const
{
    assert(_bF);
    return _bead_forces[bead];
}

inline void Configuration::HasPos(bool b)
{
    _bPos=true;
    if(_bead_positions.size()==0 && _bPos)
        _bead_positions.resize(_nBeads);
}

inline void Configuration::HasVel(bool b)
{
    _bVel=true;
    if(_bead_velocities.size()==0 && _bVel)
        _bead_velocities.resize(_nBeads);
}

inline void Configuration::HasF(bool b)
{
    _bF=true;
    if(_bead_forces.size()==0 && _bF)
        _bead_forces.resize(_nBeads);
}

inline void Configuration::HasU(bool b)
{
    _bU=true;
    if(_bead_u.size()==0 && _bU)
        _bead_u.resize(_nBeads);
}

inline void Configuration::HasV(bool b)
{
    _bV=true;
    if(_bead_v.size()==0 && _bV)
        _bead_v.resize(_nBeads);
}

inline void Configuration::HasW(bool b)
{
    _bW=true;
    if(_bead_w.size()==0 && _bW)
        _bead_w.resize(_nBeads);
}


#endif	/* _configuration_H */

