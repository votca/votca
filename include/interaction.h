/// \addtogroup csg
///@{
// 
// File:   interaction.h
// Author: ruehle
//
// Created on April 23, 2007, 2:26 PM
//

#ifndef _interaction_H
#define	_interaction_H

#include "configuration.h"

#include <string>
#include <sstream>
using namespace std;

/**
    \brief base calss for all interactions

    This is the base class for all interactions.

    \todo double names/groups right, add molecules!!
*/
class Interaction
{
public:
    virtual ~Interaction() { }
    virtual double EvaluateVar(const Configuration &conf) = 0;
        
    void setName(const string &name) { _name = name; }
    string getName() const { return _name; }
        
    void setGroup(const string &group) { _group = group; RebuildName(); }
    const string &getGroup() const { return _group; }

    // the group id is set by topology, when interaction is added to it
    // \todo if the group name is changed later, group id should be updated by topology
    int getGroupId() { return _group_id; }
    void setGroupId(int id) { _group_id = id; }

    void setIndex(const int &index) { _index = index; RebuildName(); }
    const int &getIndex() const { return _index; }
    
    void setMolecule(const int &mol) { _mol = mol; RebuildName(); }
    const int &getMolecule() const { return _mol; }
    
    virtual vec Grad(const Configuration &conf, int bead) = 0;
    int BeadCount() { return _beads.size(); }
    int getBeadId(int bead) { return _beads[bead]; }
    
protected:
    int _index;
    string _group;
    int _group_id;
    string _name;
    int _mol;
    vector<int> _beads;
    
    void RebuildName();
};

inline void Interaction::RebuildName()
{
    stringstream s;
    s << _mol+1  << ":" << _group << ":" << _index+1;
    _name = s.str();
}

/**
    \brief bond interaction
*/
class IBond : public Interaction
{
public:
    IBond(int bead1, int bead2)
        { _beads.resize(2); _beads[0] = bead1; _beads[1] = bead2; }
    
    double EvaluateVar(const Configuration &conf);
    vec Grad(const Configuration &conf, int bead);

private:
};

/**
    \brief angle interaction
*/
class IAngle : public Interaction
{
public:
    IAngle(int bead1, int bead2, int bead3)
        { _beads.resize(3); _beads[0] = bead1; _beads[1] = bead2; _beads[2] = bead3;}

    double EvaluateVar(const Configuration &conf);
    vec Grad(const Configuration &conf, int bead);

private:
};

/**
    \brief dihedral interaction
*/
class IDihedral : public Interaction
{
public:
    IDihedral(int bead1, int bead2, int bead3, int bead4)
        { _beads.resize(4); _beads[0] = bead1; _beads[1] = bead2; _beads[2] = bead3; _beads[3] = bead4;}
        
    double EvaluateVar(const Configuration &conf);
    vec Grad(const Configuration &conf, int bead) {}

private:
};

inline double IBond::EvaluateVar(const Configuration &conf)
{
    return abs(conf.getDist(_beads[0], _beads[1]));
}

inline vec IBond::Grad(const Configuration &conf, int bead)
{
    vec r = conf.getDist(_beads[0], _beads[1]); 
    r.normalize();
    return (bead == 0) ? -r : r;
}
    
inline double IAngle::EvaluateVar(const Configuration &conf)
{
    vec v1(conf.getDist(_beads[1], _beads[0]));
    vec v2(conf.getDist(_beads[1], _beads[2]));
    return  acos(v1*v2/sqrt((v1*v1) * (v2*v2)));    
}

inline vec IAngle::Grad(const Configuration &conf, int bead)
{
    vec v1(conf.getDist(_beads[1], _beads[0]));
    vec v2(conf.getDist(_beads[1], _beads[2]));
    
    double acos_prime = -1.0 / (sqrt(1 - (v1 * v2)/( abs(v1) * abs(v2) ) ));
    switch (bead) {
        case (0): return acos_prime * (-v2 / ( abs(v1)*abs(v2) ) +  (v1*v2) * v1 / ( abs(v2)*abs(v1)*abs(v1)*abs(v1) ) ); break;
        case (1): return acos_prime * ( (v1+v2)/(abs(v1) * abs(v2)) - (v1 * v2) * ((v2*v2) * v1 + (v1*v1) * v2 ) / ( abs(v1)*abs(v1)*abs(v1)*abs(v2)*abs(v2)*abs(v2) ) ); break;
        case (2): return acos_prime * (-v1 / ( abs(v1)*abs(v2) ) +  (v1*v2) * v2 / ( abs(v1)*abs(v2)*abs(v2)*abs(v2) ) ); break;
    }
}

inline double IDihedral::EvaluateVar(const Configuration &conf)
{
    vec v1(conf.getPos(_beads[1]) - conf.getPos(_beads[0]));
    vec v2(conf.getPos(_beads[2]) - conf.getPos(_beads[1]));
    vec v3(conf.getPos(_beads[3]) - conf.getPos(_beads[2]));
    vec n1, n2;
    n1 = v1^v2; // calculate the normal vector
    n2 = v2^v3; // calculate the normal vector
    double sign = (v1*n2 < 0) ? -1 : 1;
    return sign*acos(n1*n2/sqrt((n1*n1) * (n2*n2)));
    //return sign*acos(n1*n2/sqrt((n1*n1) * (n2*n2))) + 1;
    //return pow(acos(n1*n2/sqrt((n1*n1) * (n2*n2))), 2);
}

#endif	/* _interaction_H */

/// @}
