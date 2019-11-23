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
#pragma once
#ifndef VOTCA_CSG_BASEBEAD_H
#define VOTCA_CSG_BASEBEAD_H

#include <assert.h>
#include <memory>
#include <votca/csg/moleculeitem.h>
#include <votca/csg/topologyitem.h>
#include <votca/tools/constants.h>
#include <votca/tools/eigen.h>
#include <votca/tools/name.h>
#include <votca/tools/types.h>
namespace TOOLS = votca::tools;

namespace votca {
namespace csg {

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
  virtual ~BaseBead() = default;

  /// Gets the id of the bead
  int getId() const noexcept { return id_; }

  /// Sets the id of the bead
  void setId(const int &id) noexcept { id_ = id; }

  /// Gets the name of the bead
  std::string getName() const { return name_.getName(); }

  /// Sets the name of the bead
  void setName(std::string name) { return name_.setName(name); }

  /// Sets the molecule the bead is attached too
  // void setMolecule(Molecule *molecule) {
  // molecule_item_.setMolecule(molecule); }
  void setMoleculeId(const int &molecule_id) noexcept {
    molecule_id_ = molecule_id;
  }

  /// Gets the molecule pointer the bead is attached too
  // Molecule *getMolecule() const { return molecule_item_.getMolecule(); }
  int getMoleculeId() const noexcept { return molecule_id_; }

  /// Gets the topology pointer the bead is attached too
  Topology *getParent() const { return topology_item_.getParent(); }

  /**
   * get the bead type
   * \return const string
   */
  virtual const std::string getType() const noexcept { return type_; }

  /**
   * set the bead type
   * \param bead type object
   */
  virtual void setType(const std::string &type) noexcept { type_ = type; }

  /**
   * @brief Returns the element type of the bead
   *
   * @return either the element symbol i.e. "Si" for silcon or unassigned if it
   * has not been specified.
   */
  std::string getElement() const noexcept { return element_symbol_; }

  /**
   * get the mass of the base bead
   * \return - base bead mass
   */
  virtual const double &getMass() const noexcept { return mass_; }

  /**
   * set the mass of the base bead
   * \param - base bead mass
   */
  virtual void setMass(const double &m) noexcept { mass_ = m; }

  /**
   * set the position of the base bead
   * \param - base bead position
   */
  virtual void setPos(const Eigen::Vector3d &bead_position);

  /**
   * get the position of the base bead
   * \return base bead position
   */
  virtual const Eigen::Vector3d &getPos() const;

  /**
   * direct access (read/write) to the position of the base bead
   * \return reference to position
   */
  virtual Eigen::Vector3d &Pos() {
    assert(bead_position_set_ && "Position is not set.");
    return bead_position_;
  }

  /** does this configuration store positions? */
  bool HasPos() const noexcept { return bead_position_set_; }

  /** set has position to true */
  void HasPos(const bool &true_or_false) noexcept {
    bead_position_set_ = true_or_false;
  }

 protected:
  BaseBead() : topology_item_(nullptr), mass_(0.0), bead_position_set_(false){};

  TopologyItem topology_item_;

  std::string type_ = tools::topology_constants::unassigned_bead_type;
  int id_ = tools::topology_constants::unassigned_residue_id;
  int molecule_id_ = tools::topology_constants::unassigned_molecule_id;
  std::string element_symbol_ = tools::topology_constants::unassigned_element;
  TOOLS::Name name_;

  double mass_;
  Eigen::Vector3d bead_position_;

  bool bead_position_set_;
};

inline void BaseBead::setPos(const Eigen::Vector3d &bead_position) {
  bead_position_set_ = true;
  bead_position_ = bead_position;
}

inline const Eigen::Vector3d &BaseBead::getPos() const {
  assert(bead_position_set_ &&
         "Cannot get bead position as it has not been set.");
  return bead_position_;
}
}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BASEBEAD_H
