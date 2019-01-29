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

#ifndef _VOTCA_CSG_BASEBEAD_H
#define _VOTCA_CSG_BASEBEAD_H

#include <assert.h>
#include <memory>
#include <votca/csg/moleculeitem.h>
#include <votca/csg/topologyitem.h>
#include <votca/tools/identity.h>
#include <votca/tools/name.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {
using namespace votca::tools;

class BeadType;

/**
 * \brief information about a base bead
 *
 * The Base Bead class describes the core functionality of an atom or a coarse
 * grained bead. It stores information like the id, the name, the mass, the
 * charge and the residue it belongs to and the position
 *
 **/
class BaseBead {
 public:
  /**
   * destructor
   */
  virtual ~BaseBead() {}

  /// Gets the id of the bead
  int getId() const { return id_.getId(); }

  /// Sets the id of the bead
  void setId(int id) { id_.setId(id); }

  /// Gets the name of the bead
  std::string getName() const { return name_.getName(); }

  /// Sets the name of the bead
  void setName(std::string name) { return name_.setName(name); }

  /// Sets the molecule the bead is attached too
  void setMolecule(Molecule *molecule) { molecule_item_.setMolecule(molecule); }

  /// Gets the molecule pointer the bead is attached too
  Molecule *getMolecule() const { return molecule_item_.getMolecule(); }

  /// Gets the topology pointer the bead is attached too
  Topology *getParent() const { return topology_item_.getParent(); }

  /**
   * get the bead type
   * \return const bead type pointer
   */
  virtual const std::weak_ptr<BeadType> getType() const { return type_; }

  /**
   * set the bead type
   * \param bead type object
   */
  virtual void setType(std::weak_ptr<BeadType> type) { type_ = type; }

  /**
   * get the bead type
   * \return - non constant bead type pointer
   */
  virtual std::weak_ptr<BeadType> Type() const { return type_; }

  /**
   * get the name of the bead type
   * \return - string indicates the name
   **/
  std::string getBeadTypeName();

  /**
   * get the id of the bead type
   * \return - int indicated the id
   **/
  int getBeadTypeId();

  /**
   * get the mass of the base bead
   * \return - base bead mass
   */
  virtual const double &getMass() const { return mass_; }

  /**
   * set the mass of the base bead
   * \param - base bead mass
   */
  virtual void setMass(const double &m) { mass_ = m; }

  /**
   * set the position of the base bead
   * \param - base bead position
   */
  virtual void setPos(const vec &bead_position);

  /**
   * get the position of the base bead
   * \return base bead position
   */
  virtual const vec &getPos() const;

  /**
   * direct access (read/write) to the position of the base bead
   * \return reference to position
   */
  virtual vec &Pos() {
    assert(bead_position_set_ && "Position is not set.");
    return bead_position_;
  }

  /** does this configuration store positions? */
  bool HasPos() const { return bead_position_set_; }

  /** set has position to true */
  void HasPos(bool true_or_false) { bead_position_set_ = true_or_false; }

 protected:
  BaseBead()
      : topology_item_(nullptr),
        molecule_item_(nullptr),
        mass_(0.0),
        bead_position_set_(false){};

  TopologyItem topology_item_;
  MoleculeItem molecule_item_;

  Identity<int> id_;
  Name name_;
  std::weak_ptr<BeadType> type_;

  double mass_;
  vec bead_position_;

  bool bead_position_set_;
};

inline void BaseBead::setPos(const vec &bead_position) {
  bead_position_set_ = true;
  bead_position_ = bead_position;
}

inline const vec &BaseBead::getPos() const {
  assert(bead_position_set_ &&
         "Cannot get bead position as it has not been set.");
  return bead_position_;
}
}  // namespace csg
}  // namespace votca

#endif  // _VOTCA_CSG_BASEBEAD_H
