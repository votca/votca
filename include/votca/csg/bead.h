/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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
#define _VOTCA_CSG_BEAD_H

#include "basebead.h"
#include <cassert>
#include <string>
#include <votca/tools/property.h>
#include <votca/tools/types.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {

class BeadType;
class Topology;
class Molecule;

using namespace votca::tools;

/**
 * \brief information about a bead
 *
 * The Bead class describes an atom or a coarse grained bead. It stores
 * information like the id, the name, the mass, the
 * charge and the residue it belongs to. The coordinates are stored in the
 * configuration class.
 *
 * \todo change resnr to pointer
 * \todo make sure bead belongs to topology
 **/
class Bead : public BaseBead {
 public:
  /**
   * destructor
   */
  virtual ~Bead() {}

  /**
   * get the residu number of the bead
   * \return residue id
   */
  const int &getResnr() const { return residue_number_; }

  /**
   * get the mass of the bead
   * \return bead mass
   */
  const double &getM() const {
    std::cerr << "WARNING getM is depricated use getMass" << std::endl;
    return getMass();
  }

  /**
   * set the mass of the bead
   * \param m bead mass
   */
  void setM(const double &m) {
    std::cerr << "WARNING setM is depricated use setMass" << std::endl;
    setMass(m);
  }

  /**
   * get the charge of the bead
   * \return - base bead charge
   */
  virtual const double &getQ() const { return charge_; }

  /**
   * set the charge of the base bead
   * \param[in] - base bead position
   */
  virtual void setQ(const double &q) { charge_ = q; }

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
  byte_t getSymmetry() const { return symmetry_; }

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
   * axes are denoted as u, v and w. Currently the non-spherical mapping is
   * hardcoded and
   * the axis u is calculated by the eigenvector with the lowest eigenvector of
   * the mapped beads and has the meaning of a normal vector if the reference
   * beads
   * have a disc like shape. The sign of the normal vector is determined in
   * combination
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
   * axes are denoted as u, v and w. Currently the non-spherical mapping is
   * hardcoded and
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
   * axes are denoted as u, v and w. Currently the non-spherical mapping is
   * hardcoded and
   * the axis w is orthogonal to u and v.
   *
   * \return bead orientation w
   */
  const vec &getW() const;

  /**
   * direct access (read/write) to the velocity of the bead
   * \return reference to velocity
   */
  vec &Vel() {
    assert(bead_velocity_set_);
    return velocity_;
  }

  /**
   * direct access (read/write) to orientation u of the bead
   * \return reference to u
   */
  vec &U() {
    assert(bU_);
    return u_;
  }

  /**
   * direct access (read/write) to the orientation v of the bead
   * \return reference to v
   */
  vec &V() {
    assert(bV_);
    return v_;
  }

  /**
   * direct access (read/write) to the orientation w of the bead
   * \return reference to w
   */
  vec &W() {
    assert(bW_);
    return w_;
  }

  /**
   * direct access (read/write) to the force of the bead
   * \return reference to force
   */
  vec &F() {
    assert(bead_force_set_);
    return bead_force_;
  }

  /**
   * set force acting on bead
   * @param bead_force force
   */
  void setF(const vec &bead_force);

  /**
   * \brief get the force acting on the bead
   *
   * Forces have to be provided by the trajectory. If beads are mapped, forces
   * of coarse-grained beads are also calculated.
   *
   * \return force on bead
   */
  const vec &getF() const;

  /** does this configuration store velocities? */
  bool HasVel() { return bead_velocity_set_; }

  /** does this configuration store forces? */
  bool HasF() { return bead_force_set_; }

  /** does this configuration store u-orientations? */
  bool HasU() { return bU_; }

  /** does this configuration store v-orientations? */
  bool HasV() { return bV_; }

  /** does this configuration store w-orientations? */
  bool HasW() { return bW_; }

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
   * If it is a mapped beads, returns te bead id the cg bead was created from
   * \return vector of bead ids of reference atoms
   */
  std::vector<int> &ParentBeads() { return parent_beads_; };

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
  template <typename T>
  void setUserData(T *userdata) {
    _userdata = (void *)userdata;
  }

  /**
   * get userdata attached to bead
   * @return pointer to userdata
   */
  template <typename T>
  T *getUserData() {
    return (T *)_userdata;
  }

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
  Property &Options() { return *options_; }

  /**
   * update pointer to options object of bead
   * \param options pointer to options object of bead
   */
  void setOptions(Property &options) { options_ = &options; }

 protected:
  std::vector<int> parent_beads_;

  // TODO: this is so far a pointer. this should change! each bead should have
  // own options.
  Property *options_;

  byte_t symmetry_;
  double charge_;

  int residue_number_;

  vec velocity_, bead_force_, u_, v_, w_;

  bool bead_velocity_set_;
  bool bU_;
  bool bV_;
  bool bW_;
  bool bead_force_set_;

  /// constructur
  Bead(Topology *owner, int id, std::weak_ptr<BeadType> type, byte_t symmetry,
       std::string name, int resnr, double m, double q)
      : symmetry_(symmetry), charge_(q), residue_number_(resnr) {
    topology_item_._parent = owner;
    setId(id);
    setType(type);
    setName(name);
    setMass(m);
    bPos_ = false;
    bead_velocity_set_ = false;
    bU_ = false;
    bV_ = false;
    bW_ = false;
    bead_force_set_ = false;
  }

  void *_userdata;

  friend class Topology;
  friend class Molecule;
};

inline void Bead::setVel(const vec &r) {
  bead_velocity_set_ = true;
  velocity_ = r;
}

inline const vec &Bead::getVel() const {
  assert(bead_velocity_set_);
  return velocity_;
}

inline void Bead::setU(const vec &u) {
  bU_ = true;
  u_ = u;
}

inline const vec &Bead::getU() const {
  assert(bU_);
  return u_;
}

inline void Bead::setV(const vec &v) {
  bV_ = true;
  v_ = v;
}

inline const vec &Bead::getV() const {
  assert(bV_);
  return v_;
}

inline void Bead::setW(const vec &w) {
  bW_ = true;
  w_ = w;
}

inline const vec &Bead::getW() const {
  assert(bW_);
  return w_;
}

inline void Bead::setF(const vec &bead_force) {
  bead_force_set_ = true;
  bead_force_ = bead_force;
}

inline const vec &Bead::getF() const {
  assert(bead_force_set_);
  return bead_force_;
}

inline void Bead::HasVel(bool b) { bead_velocity_set_ = b; }

inline void Bead::HasF(bool b) { bead_force_set_ = b; }

inline void Bead::HasU(bool b) { bU_ = b; }

inline void Bead::HasV(bool b) { bV_ = b; }

inline void Bead::HasW(bool b) { bW_ = b; }
}  // namespace csg
}  // namespace votca

#endif  // _VOTCA_CSG_BEAD_H
