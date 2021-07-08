/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_CSG_TOPOLOGY_H
#define VOTCA_CSG_TOPOLOGY_H
#pragma once

// Standard includes
#include <cassert>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

// Third party includes
#include <boost/container/deque.hpp>

// VOTCA includes
#include <votca/tools/types.h>

// Local VOTCA includes
#include "bead.h"
#include "boundarycondition.h"
#include "exclusionlist.h"
#include "molecule.h"
#include "openbox.h"
#include "orthorhombicbox.h"
#include "residue.h"
#include "triclinicbox.h"

namespace votca {
namespace csg {

class Interaction;

/* Boost deque has been chosen and contents have been replaced with objects
 * as opposed to heap allocated types:
 * 1. To get rid of indirection
 * 2. Clarify ownership
 * 3. To ensure pointers are not invalidated if the container changes size
 * that is not a guarantee with a vector
 * 4. To provide better contiguous memory access, not possible with std::deque
 * or list
 */
typedef boost::container::deque_options<
    boost::container::block_size<sizeof(Residue) * 4>>::type block_residue_x4_t;
typedef boost::container::deque_options<
    boost::container::block_size<sizeof(Molecule) * 4>>::type
    block_molecule_4x_t;
typedef boost::container::deque_options<
    boost::container::block_size<sizeof(Bead) * 4>>::type block_bead_x4_t;

using MoleculeContainer =
    boost::container::deque<Molecule, void, block_molecule_4x_t>;
using BeadContainer = boost::container::deque<Bead, void, block_bead_x4_t>;
using ResidueContainer =
    boost::container::deque<Residue, void, block_residue_x4_t>;
using InteractionContainer = std::vector<Interaction *>;

/**
 * \brief topology of the whole system
 *
 * The Topology class stores the topology of the system like the beads, bonds,
 * molecules and residues.
 *
 **/
class Topology {
 public:
  /// constructor
  Topology() { bc_ = std::make_unique<OpenBox>(); }
  ~Topology();

  /**
   * \brief Cleans up all the stored data
   */
  void Cleanup();

  /**
   * \brief Creates a new Bead
   *
   * \param[in] symmetry symmetry of the bead, 1: spherical 3: ellipsoidal
   * \param[in] name name of the bead
   * \param[in] type bead type
   * \param[in] resnr residue number
   * \param[in] m mass
   * \param[in] q charge
   * \return pointer to created bead
   *
   * The function creates a new bead and adds it to the list of beads.
   */
  Bead *CreateBead(Bead::Symmetry symmetry, std::string name, std::string type,
                   Index resnr, double m, double q);

  /**
   * \brief Creates a new molecule
   * \param[in] name name of the molecule
   * \return pointer to created molecule
   */
  Molecule *CreateMolecule(std::string name);

  /**
   *  \brief checks weather molecules with the same name really contain the same
   * number of beads
   */
  void CheckMoleculeNaming(void);

  /**
   * \brief Create a new resiude.
   *
   * @param[in] name residue name
   * @return created residue
   */
  Residue &CreateResidue(std::string name);
  Residue &CreateResidue(std::string name, Index id);

  /**
   * \brief Create molecules based on the residue.
   *
   * This function scans the topology and creates molecules based on the resiude
   * id. All beads with the same resid are put Index one molecule.
   */
  void CreateMoleculesByResidue();

  /**
   * \brief put the whole topology in one molecule
   * \param name name of the new molecule
   *
   *  This function creates one big molecule for all beads in the topology.
   */
  void CreateOneBigMolecule(std::string name);

  /**
   * \brief create molecules based on blocks of atoms
   * \param[in] name molecule name
   * \param[in] first first bead
   * \param[in] nbeads number of beads per molecule
   * \param[in] nmolecules number of molecules
   */
  void CreateMoleculesByRange(std::string name, Index first, Index nbeads,
                              Index nmolecules);

  /**
   * \brief number of molecules in the system
   * @return number of molecule in topology
   */
  Index MoleculeCount() const { return molecules_.size(); }

  /**
   * number of beads in the system
   * @return number of beads in the system
   */
  Index BeadCount() const { return beads_.size(); }

  /**
   * number of residues in the system
   * \return number of residues
   */
  Index ResidueCount() const { return residues_.size(); }

  /**
   * get molecule by index
   * @param index molecule number
   * @return pointer to molecule
   */
  Molecule *MoleculeByIndex(Index index);

  /**
   * access containter with all beads
   * @return bead container
   */
  BeadContainer &Beads() { return beads_; }

  /**
   * access containter with all residues
   * @return bead container
   */
  ResidueContainer &Residues() { return residues_; }
  const ResidueContainer &Residues() const { return residues_; }

  /**
   * access  containter with all molecules
   * @return molecule container
   */
  MoleculeContainer &Molecules() { return molecules_; }
  const MoleculeContainer &Molecules() const { return molecules_; }

  /**
   * access containter with all bonded interactions
   * @return bonded interaction container
   */
  InteractionContainer &BondedInteractions() { return interactions_; }
  const InteractionContainer &BondedInteractions() const {
    return interactions_;
  }

  void AddBondedInteraction(Interaction *ic);
  std::vector<Interaction *> InteractionsInGroup(const std::string &group);

  /**
   * \brief Determine if a bead type exists.
   *
   * @return bool true if it has been registered
   **/
  bool BeadTypeExist(std::string type) const;

  /**
   * \brief Register the bead type with the topology object.
   *
   * Records are kept of the different bead types in the topology object. This
   * method stores the bead type.
   **/
  void RegisterBeadType(std::string type);

  /**
   * \brief Given a bead type this method returns the id associated with the
   * type
   *
   * @param[in] string name of the type
   * @return Index the id of the type
   **/
  Index getBeadTypeId(std::string type) const;

  /**
   * \brief Returns a pointer to the bead with index i
   *
   * @param[in] Index i is the id of the bead
   * @return Bead * is a pointer to the bead
   **/
  Bead *getBead(const Index i) { return &beads_[i]; }
  const Bead *getBead(const Index i) const { return &beads_[i]; }
  Residue &getResidue(const Index i) { return residues_[i]; }
  const Residue &getResidue(const Index i) const { return residues_[i]; }
  Molecule *getMolecule(const Index i) { return &molecules_[i]; }
  const Molecule *getMolecule(const Index i) const { return &molecules_[i]; }

  /**
   * delete all molecule information
   */
  void ClearMoleculeList() { molecules_.clear(); }

  /**
   * \brief adds all the beads+molecules+residues from other topology
   * \param top topology to add
   */
  void Add(Topology *top);

  /**
   * \brief copy topology data of different topology
   * \param top topology to copy from
   */
  void CopyTopologyData(Topology *top);

  /**
   *  \brief rename all the molecules in range
   * \param range range string of type 1:2:10 = 1, 3, 5, 7, ...
   * \param name new name of molecule
   * range is a string which is parsed by RangeParser,
   */
  void RenameMolecules(std::string range, std::string name);

  /**
   *  \brief rename all the bead types
   * \param name current rame of the bead type
   * \param newname new name of bead type
   */
  void RenameBeadType(std::string name, std::string newname);

  /**
   *  \brief set the mass of all the beads of a certain type
   * \param name the bead type
   * \param value mass value
   */
  void SetBeadTypeMass(std::string name, double value);

  /**
   * set the simulation box
   * \param box triclinic box matrix
   */
  void setBox(const Eigen::Matrix3d &box, BoundaryCondition::eBoxtype boxtype =
                                              BoundaryCondition::typeAuto) {
    // determine box type automatically in case boxtype==typeAuto
    if (boxtype == BoundaryCondition::typeAuto) {
      boxtype = autoDetectBoxType(box);
    }

    switch (boxtype) {
      case BoundaryCondition::typeTriclinic:
        bc_ = std::make_unique<TriclinicBox>();
        break;
      case BoundaryCondition::typeOrthorhombic:
        bc_ = std::make_unique<OrthorhombicBox>();
        break;
      default:
        bc_ = std::make_unique<OpenBox>();
        break;
    }

    bc_->setBox(box);
  };

  /**
   * get the simulation box
   * \return triclinic box matrix
   */
  const Eigen::Matrix3d &getBox() const { return bc_->getBox(); };

  /**
   * @brief Return the boundary condition object
   */
  const BoundaryCondition &getBoundary() const {
    assert(bc_ && "Cannot return boundary condition is null");
    return *bc_;
  };
  /**
   * set the time of current frame
   * \param t simulation time in ns
   */
  void setTime(double t) { time_ = t; };

  /**
   * get the time of current frame
   * \return simulation time in ns
   */
  double getTime() const { return time_; };

  /**
   * set the step number of current frame
   * \param s step number
   */
  void setStep(Index s) { step_ = s; };

  /**
   * get the step number of current frame
   * \return step number
   */
  Index getStep() const { return step_; };

  /**
   * Sets the particle group. (For the H5MD file format)
   * \param particle_group The name of a particle group.
   */
  void setParticleGroup(std::string particle_group) {
    particle_group_ = particle_group;
  };

  /**
   * Gets the particle group.
   * \return The name of a particle group.
   */
  std::string getParticleGroup() const { return particle_group_; };

  /**
   * \brief pbc correct distance of two beads
   * \param bead1 index of first bead
   * \param bead2 index of second bead
   * \return distance vector
   *
   * calculates the smallest distance between two beads with correct treatment
   * of pbc
   */
  Eigen::Vector3d getDist(Index bead1, Index bead2) const;

  /**
   * \brief calculate shortest vector connecting two points
   * \param r1 first point
   * \param r2 second point
   * \return distance vector
   *
   * calculates the smallest distance between two points with correct treatment
   * of pbc
   */
  Eigen::Vector3d BCShortestConnection(const Eigen::Vector3d &r_i,
                                       const Eigen::Vector3d &r_j) const;

  /**
   * \brief return the shortest box size
   * \return shortest size
   *
   * Calculates the shortest length to connect two sides of the box
   */
  double ShortestBoxSize() const;

  /**
   *  calculates the box volume
   *  \return box volume
   */
  double BoxVolume() const;

  /**
   *  rebuild exclusion list
   */
  void RebuildExclusions();

  /**
   * access exclusion list
   * \return exclusion list
   */
  ExclusionList &getExclusions() { return exclusions_; }
  const ExclusionList &getExclusions() const { return exclusions_; }

  BoundaryCondition::eBoxtype getBoxType() const { return bc_->getBoxType(); }

  template <typename iteratable>
  void InsertExclusion(Bead *bead1, iteratable &l);

  bool HasVel() { return has_vel_; }
  void SetHasVel(const bool v) { has_vel_ = v; }

  bool HasForce() { return has_force_; }
  void SetHasForce(const bool v) { has_force_ = v; }

 protected:
  std::unique_ptr<BoundaryCondition> bc_;

  BoundaryCondition::eBoxtype autoDetectBoxType(
      const Eigen::Matrix3d &box) const;

  /// bead types in the topology
  std::unordered_map<std::string, Index> beadtypes_;

  /// beads in the topology
  BeadContainer beads_;

  /// molecules in the topology
  MoleculeContainer molecules_;

  /// residues in the topology
  ResidueContainer residues_;

  /// bonded interactions in the topology
  InteractionContainer interactions_;

  ExclusionList exclusions_;

  std::map<std::string, Index> interaction_groups_;

  std::map<std::string, std::vector<Interaction *>> interactions_by_group_;

  double time_ = 0.0;
  Index step_ = 0;
  bool has_vel_ = false;
  bool has_force_ = false;

  /// The particle group (For H5MD file format)
  std::string particle_group_ = "unassigned";
};

inline Bead *Topology::CreateBead(Bead::Symmetry symmetry, std::string name,
                                  std::string type, Index resnr, double m,
                                  double q) {

  beads_.push_back(Bead(beads_.size(), type, symmetry, name, resnr, m, q));
  return &beads_.back();
}

inline Molecule *Topology::CreateMolecule(std::string name) {
  molecules_.push_back(Molecule(molecules_.size(), name));
  return &molecules_.back();
}

inline Residue &Topology::CreateResidue(std::string name, Index id) {
  // Note that Residue constructor is intentionally private and only topology
  // class can create it, hence emplace back will not work because the vector
  // class does not have access to the constructor.
  residues_.push_back(Residue(id, name));
  return residues_.back();
}

inline Residue &Topology::CreateResidue(std::string name) {
  // Note that Residue constructor is intentionally private and only topology
  // class can create it, hence emplace back will not work because the vector
  // class does not have access to the constructor.
  residues_.push_back(Residue(residues_.size(), name));
  return residues_.back();
}

inline Molecule *Topology::MoleculeByIndex(Index index) {
  return &molecules_[index];
}

template <typename iteratable>
inline void Topology::InsertExclusion(Bead *bead1, iteratable &l) {
  exclusions_.InsertExclusion(bead1, l);
}

}  // namespace csg
}  // namespace votca

#include "interaction.h"

#endif  // VOTCA_CSG_TOPOLOGY_H
