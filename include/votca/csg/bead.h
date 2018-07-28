/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_CSG_BEAD_H
#define	_VOTCA_CSG_BEAD_H

#include <string>
#include <votca/tools/types.h>
#include <votca/tools/vec.h>
#include <votca/tools/property.h>
#include <assert.h>
#include "beadtype.h"
#include "topologyitem.h"

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;
class Molecule;

/**
    \brief information about a bead
 
    The Bead class describes an atom or a coarse grained bead. It stores information like the id, the name, the mass, the
    charge and the residue it belongs to. The coordinates are stored in the configuration class.

    \todo change resnr to pointer
    \todo make sure bead belongs to topology
*/
class Bead : public TopologyItem
{
public:   

    /**
     * destructor
     */
    virtual ~Bead() {}

    /**
     * get the id of the bead
     *
     * \return bead id
     */
    const int &getId() const { return _id; }
    
    /**
     * get bead name
     * \return bead name
     */
    const string &getName() const { return _name; }
    
    /**
     * set bead name
     * \param name bead name
     */
    void setName(const string &name) { _name=name; }

    /**
     * get the bead type
     * \return bead type object
     */
    const BeadType *getType() const { return _type; }

    /**
     * set the bead type
     * \param bead type object
     */
    void setType(BeadType *type) { _type=type; }

    /**
     * get the bead type pointer (not constant)
     * \return bead type object
     */
    BeadType *Type() const { return _type; }

    /**
     * get the residu number of the bead
     * \return residue id
     */
    const int &getResnr() const { return _resnr; }

    /**
     * get the mass of the bead
     * \return bead mass
     */
    const double &getM() const { return _m; }

    /**
     * get the charge of the bead
     * \return bead charge
     */
    const double &getQ() const { return _q; }
    
    /**
     * set the mass of the bead
     * \param m bead mass
     */
    void setM(const double &m) { _m=m; }

    /**
     * set the charge of the bead
     * \param q bead charge
     */
    void setQ(const double &q) { _q=q; }

    /**
     * \brief get the symmetry of the bead
     *
     * Returns the number of unique axis of the bead, it can be
     * 1 for a spherical bead
     * 3 for an ellipsoidal bead
     * 2 (currently not used), could be disk like particle
     *
     * \return bead symmetry
     */
    byte_t getSymmetry() const { return _symmetry; }

    /**
     * set the position of the bead
     * \param r bead position
     */
    void setPos(const vec &r);

    /**
     * get the position of the bead
     * \return bead position
     */
    const vec &getPos() const;

    /**
     * set the velocity of the bead
     * @param r bead velocity
     */
    void setVel(const vec &r);

    /**
     * get the velocity of the bead
     * \return bead velocity
     */
    const vec &getVel() const;
    
    /**
     * \brief set first orientation (normal vector) vector of bead
     *
     * see getU for details
     *
     * @param u bead orientation u
     */
    void setU(const vec &u);

    /**
     * \brief get first orientation (normal vector) vector of bead
     *
     * Non-spherical beads (symmetry 3) have a internal coordinates system and the
     * axes are denoted as u, v and w. Currently the non-spherical mapping is hardcoded and
     * the axis u is calculated by the eigenvector with the lowest eigenvector of
     * the mapped beads and has the meaning of a normal vector if the reference beads
     * have a disc like shape. The sign of the normal vector is determined in combination
     * with the vectors v and w to build up a right handed (??) coordinate system.
     *
     * \return bead orientation u
     */
    const vec &getU() const;

    /**
     * \brief set second orientation vector of bead
     *
     * see getV for details
     *
     * @param v bead orientation v
     */
    void setV(const vec &v);

    /**
     * \brief get second orientation vector of bead
     *
     * Non-spherical beads (symmetry 3) have a internal coordinates system and the
     * axes are denoted as u, v and w. Currently the non-spherical mapping is hardcoded and
     * the axis v is the vector which connects first and second reference atom
     * in the mapping (only orthogonal component to u).
     *
     * \return bead orientation u
     */
    const vec &getV() const;

    /**
     * \brief set third orientation vector of bead
     *
     * see getW for details
     *
     * @param w bead orientation w
     */
    void setW(const vec &w);

    /**
     * \brief get third orientation vector of bead
     *
     * Non-spherical beads (symmetry 3) have a internal coordinates system and the
     * axes are denoted as u, v and w. Currently the non-spherical mapping is hardcoded and
     * the axis w is orthogonal to u and v.
     *
     * \return bead orientation w
     */
    const vec &getW() const;
        
    /**
     * direct access (read/write) to the position of the bead
     * \return reference to position 
     */
    vec &Pos() { return _pos; }

    /**
     * direct access (read/write) to the velocity of the bead
     * \return reference to velocity
     */
    vec &Vel() { return _vel; }

    /**
     * direct access (read/write) to orientation u of the bead
     * \return reference to u
     */
    vec &U() { return _u; }

    /**
     * direct access (read/write) to the orientation v of the bead
     * \return reference to v
     */
    vec &V() { return _v; }

    /**
     * direct access (read/write) to the orientation w of the bead
     * \return reference to w
     */
    vec &W() { return _w; }

    /**
     * direct access (read/write) to the force of the bead
     * \return reference to force
     */
    vec &F() { return _f; }
 
    /**
     * set force acting on bead
     * @param F force
     */
    void setF(const vec &F);

    /**
     * \brief get the force acting on the bead
     *
     * Forces have to be provided by the trajectory. If beads are mapped, forces
     * of coarse-grained beads are also calculated.
     *
     * \return force on bead
     */
    const vec &getF() const;

    /** does this configuration store positions? */
    bool HasPos() {return _bPos; }
    
    /** does this configuration store velocities? */
    bool HasVel() {return _bVel; }
    
    /** does this configuration store forces? */
    bool HasF() {return _bF; }
    
    /** does this configuration store u-orientations? */
    bool HasU() {return _bU; }
    
    /** does this configuration store v-orientations? */
    bool HasV() {return _bV; }
    
    /** does this configuration store w-orientations? */
    bool HasW() {return _bW; }
        
    /** dos the bead store a position */
    void HasPos(bool b);

    /** dos the bead store a velocity */
    void HasVel(bool b);
    
    /** dos the bead store a force */
    void HasF(bool b);

    /** doe the bead store an orientation u */
    void HasU(bool b);

    /** doe the bead store an orientation v */
    void HasV(bool b);

    /** doe the bead store an orientation w */
    void HasW(bool b);
    
    /**
     * molecule the bead belongs to
     * \return Molecule object
     */
    Molecule *getMolecule() { return _mol; }

    void setMolecule( Molecule * mol);
    /**
     * If it is a mapped beads, returns te bead id the cg bead was created from
     * \return vector of bead ids of reference atoms
     */
    vector<int> &ParentBeads() { return _parent_beads; };

    /**
     * \brief Function to add arbitrary user data to bead
     *
     * The user can attach pointers to own created objects to beads. Currently
     * the user has to take care about deletion of the objects at the end.
     *
     * \todo change this to shared_pointer
     *
     * \param userdata userdata
     */
    template<typename T>
    void setUserData(T *userdata) { _userdata = (void*)userdata; }

    /**
     * get userdata attached to bead
     * @return pointer to userdata
     */
    template<typename T>
    T *getUserData() { return (T *)_userdata; }

    /**
     * \brief Additional options of bead
     *
     * The options object stores additional options which can be attached to
     * the bead. For mapped beads, it contains all the values which were specified
     * in the xml mapping file. This allows to at user defined options to the xml
     * which are automatically read in on creation of the coare-grained bead.
     *
     * \return Property object containing options
     */
    Property &Options() { return *_options; }

    /**
     * update pointer to options object of bead
     * \param options pointer to options object of bead
     */
    void setOptions(Property &options) { _options = &options; }

protected:
    int _id;
    vector<int> _parent_beads;
    BeadType *_type;
    Molecule *_mol;
    
    // TODO: this is so far a pointer. this should change! each bead should have own options.
    Property *_options;

    byte_t _symmetry;
    string _name;
    
    int _resnr;
    
    double _m;
    double _q;
    
    vec _pos, _vel, _f, _u, _v, _w;
    
    bool _bPos;
    bool _bVel;
    bool _bU;
    bool _bV;
    bool _bW;
    bool _bF;
    
    /// constructur
    Bead(Topology *owner, int id, BeadType *type, byte_t symmetry, string name, int resnr, double m, double q)
        : TopologyItem(owner), _id(id), _type(type), _symmetry(symmetry), _name(name), _resnr(resnr), _m(m), _q(q)
    {_bPos=false;
    _bVel=false;
    _bU=false;
    _bV=false;
    _bW=false;
    _bF=false;}

    void *_userdata;
    
    friend class Topology;

    friend class Molecule;
};

inline void Bead::setMolecule(Molecule *mol){
	_mol = mol;
}

inline void Bead::setPos(const vec &r)
{
    _bPos=true;
    _pos = r;
}

inline const vec &Bead::getPos() const
{
    assert(_bPos);
    return _pos;
}

inline void Bead::setVel(const vec &r)
{
   _bVel=true;
   _vel = r;
}

inline const vec &Bead::getVel() const
{
    assert(_bVel);
    return _vel;
}

inline void Bead::setU(const vec &u)
{
    _bU=true;
    _u = u;
}

inline const vec &Bead::getU() const
{
    assert(_bU);
    return _u;
}

inline void Bead::setV(const vec &v)
{
    _bV=true;
    _v = v;
}

inline const vec &Bead::getV() const
{
    assert(_bV);
    return _v;
}

inline void Bead::setW(const vec &w)
{
    _bW=true;
    _w = w;
}

inline const vec &Bead::getW() const
{
    assert(_bW);
    return _w;
}

inline void Bead::setF(const vec &F)
{
    _bF=true;
    _f = F;
}

inline const vec &Bead::getF() const
{
    assert(_bF);
    return _f;
}

inline void Bead::HasPos(bool b)
{
    _bPos=b;
}

inline void Bead::HasVel(bool b)
{
    _bVel=b;
}

inline void Bead::HasF(bool b)
{
    _bF=b;
}

inline void Bead::HasU(bool b)
{
    _bU=b;
}

inline void Bead::HasV(bool b)
{
    _bV=b;
}

inline void Bead::HasW(bool b)
{
    _bW=b;
}

}}

#endif	// _VOTCA_CSG_BEAD_H

