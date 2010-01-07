/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _interaction_H
#define	_interaction_H

#include <string>
#include <sstream>
using namespace std;

#include "topology.h"
#include "bead.h"

/**
    \brief base calss for all interactions

    This is the base class for all interactions.

    \todo double names/groups right, add molecules!!
*/
class Interaction
{
public:
    virtual ~Interaction() { }
    virtual double EvaluateVar(const Topology &top) = 0;
        
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
    
    virtual vec Grad(const Topology &top, int bead) = 0;
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

    IBond(list<int> &beads)
        { _beads.resize(2); for(int i=0; i<2; ++i) { _beads[i] = beads.front(); beads.pop_front(); }}
    double EvaluateVar(const Topology &top);
    vec Grad(const Topology &top, int bead);

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
    IAngle(list<int> &beads)
        { _beads.resize(3); for(int i=0; i<3; ++i) { _beads[i] = beads.front(); beads.pop_front(); }}

    double EvaluateVar(const Topology &top);
    vec Grad(const Topology &top, int bead);

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
    IDihedral(list<int> &beads)
        { _beads.resize(4); for(int i=0; i<4; ++i) { _beads[i] = beads.front(); beads.pop_front(); }}
   
    double EvaluateVar(const Topology &top);
    vec Grad(const Topology &top, int bead) {}

private:
};

inline double IBond::EvaluateVar(const Topology &top)
{
    return abs(top.getDist(_beads[0], _beads[1]));
}

inline vec IBond::Grad(const Topology &top, int bead)
{
    vec r = top.getDist(_beads[0], _beads[1]); 
    r.normalize();
    return (bead == 0) ? -r : r;
}
    
inline double IAngle::EvaluateVar(const Topology &top)
{
    vec v1(top.getDist(_beads[1], _beads[0]));
    vec v2(top.getDist(_beads[1], _beads[2]));
    return  acos(v1*v2/sqrt((v1*v1) * (v2*v2)));    
}

inline vec IAngle::Grad(const Topology &top, int bead)
{
    /*vec v1(top.getDist(_beads[1], _beads[0]));
    vec v2(top.getDist(_beads[1], _beads[2]));
    
    double av1 = abs(v1);
    double av2 = abs(v2);
    double cosphi = v1*v2 / (av1*av2);
    double acos_prime = -1.0 / (sqrt(1 - (cosphi*cosphi) ));
    
    switch (bead) {
        case (0): return acos_prime * (v2 / (av1*av2) - v1*cosphi/(av1*av1)); break;
        case (1): return -acos_prime * ( (v1+v2)/(av1 * av2)) - 
                cosphi * ( v1/(av1*av1) + v2/(av2*av2) ); break;
        case (2): return acos_prime * (v1 / (av1*av2) - v2*cosphi/(av2*av2)); break;
    }
    return 0;
    /**/
    vec v1(top.getDist(_beads[1], _beads[0]));
    vec v2(top.getDist(_beads[1], _beads[2]));
    
    double acos_prime = 1.0 / (sqrt(1 - (v1*v2) * (v1*v2)/( abs(v1) * abs(v2) * abs(v1) * abs(v2) ) ));
    switch (bead) {
        case (0): return acos_prime * (-v2 / ( abs(v1)*abs(v2) ) +  (v1*v2) * v1 / ( abs(v2)*abs(v1)*abs(v1)*abs(v1) ) ); break;
        case (1): return acos_prime * ( (v1+v2)/(abs(v1) * abs(v2)) - (v1 * v2) * ((v2*v2) * v1 + (v1*v1) * v2 ) / ( abs(v1)*abs(v1)*abs(v1)*abs(v2)*abs(v2)*abs(v2) ) ); break;
        case (2): return acos_prime * (-v1 / ( abs(v1)*abs(v2) ) +  (v1*v2) * v2 / ( abs(v1)*abs(v2)*abs(v2)*abs(v2) ) ); break;
    }/**/
}

inline double IDihedral::EvaluateVar(const Topology &top)
{
    vec v1(top.getDist(_beads[0], _beads[1]));
    vec v2(top.getDist(_beads[1], _beads[2]));
    vec v3(top.getDist(_beads[2], _beads[3]));
    vec n1, n2;
    n1 = v1^v2; // calculate the normal vector
    n2 = v2^v3; // calculate the normal vector
    double sign = (v1*n2 < 0) ? -1 : 1;
    return sign*acos(n1*n2/sqrt((n1*n1) * (n2*n2)));
    //return sign*acos(n1*n2/sqrt((n1*n1) * (n2*n2))) + 1;
    //return pow(acos(n1*n2/sqrt((n1*n1) * (n2*n2))), 2);
}

#endif	/* _interaction_H */

