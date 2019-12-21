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
#ifndef VOTCA_CSG_TOPOLOGY_H
#define VOTCA_CSG_TOPOLOGY_H

#include <cassert>
#include <list>
#include <map>
#include <unordered_map>
#include <vector>

#include "bead.h"
#include "boundarycondition.h"
#include "exclusionlist.h"
#include "molecule.h"
#include "openbox.h"
#include "orthorhombicbox.h"
#include "residue.h"
#include "triclinicbox.h"
#include <votca/tools/types.h>

namespace votca {
namespace csg {

class Interaction;

using MoleculeContainer = std::vector<Molecule *>;
using BeadContainer = std::vector<Bead *>;
using ResidueContainer = std::vector<Residue *>;
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
  Topology() { _bc = new OpenBox(); }
  virtual ~Topology();

  /**
   * \brief Cleans up all the stored data
   */
  virtual void Cleanup();

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
  virtual Bead *CreateBead(Bead::Symmetry symmetry, std::string name,
                           std::string type, Index resnr, double m, double q);

  /**
   * \brief Creates a new molecule
   * \param[in] name name of the molecule
   * \return pointer to created molecule
   */
  virtual Molecule *CreateMolecule(std::string name);

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
  virtual Residue *CreateResidue(std::string name);
  virtual Residue *CreateResidue(std::string name, Index id);

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
  Index MoleculeCount() const { return _molecules.size(); }

  /**
   * number of beads in the system
   * @return number of beads in the system
   */
  Index BeadCount() const { return _beads.size(); }

  /**
   * number of residues in the system
   * \return number of residues
   */
  Index ResidueCount() const { return _residues.size(); }

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
  BeadContainer &Beads() { return _beads; }

  /**
   * access containter with all residues
   * @return bead container
   */
  ResidueContainer &Residues() { return _residues; }

  /**
   * access  containter with all molecules
   * @return molecule container
   */
  MoleculeContainer &Molecules() { return _molecules; }
  const MoleculeContainer &Molecules() const { return _molecules; }

  /**
   * access containter with all bonded interactions
   * @return bonded interaction container
   */
  InteractionContainer &BondedInteractions() { return _interactions; }
  const InteractionContainer &BondedInteractions() const { return _interactions; }

  void AddBondedInteraction(Interaction *ic);
  std::list<Interaction *> InteractionsInGroup(const std::string &group);

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
  void RegisterBeadType(std::string name);

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
  Bead *getBead(const Index i) const { return _beads[i]; }
  Residue *getResidue(const Index i) const { return _residues[i]; }
  Molecule *getMolecule(const Index i) const { return _molecules[i]; }

  /**
   * delete all molecule information
   */
  void ClearMoleculeList() { _molecules.clear(); }

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

    if (_bc) {
      delete (_bc);
    }

    switch (boxtype) {
      case BoundaryCondition::typeTriclinic:
        _bc = new TriclinicBox();
        break;
      case BoundaryCondition::typeOrthorhombic:
        _bc = new OrthorhombicBox();
        break;
      default:
        _bc = new OpenBox();
        break;
    }

    _bc->setBox(box);
  };

  /**
   * get the simulation box
   * \return triclinic box matrix
   */
  const Eigen::Matrix3d &getBox() const { return _bc->getBox(); };

  /**
   * set the time of current frame
   * \param t simulation time in ns
   */
  void setTime(double t) { _time = t; };

  /**
   * get the time of current frame
   * \return simulation time in ns
   */
  double getTime() const { return _time; };

  /**
   * set the step number of current frame
   * \param s step number
   */
  void setStep(Index s) { _step = s; };

  /**
   * get the step number of current frame
   * \return step number
   */
  Index getStep() const { return _step; };

  /**
   * Sets the particle group. (For the H5MD file format)
   * \param particle_group The name of a particle group.
   */
  void setParticleGroup(std::string particle_group) {
    _particle_group = particle_group;
  };

  /**
   * Gets the particle group.
   * \return The name of a particle group.
   */
  std::string getParticleGroup() const { return _particle_group; };

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
  Eigen::Vector3d BCShortestConnection(const Eigen::Vector3d &r1,
                                       const Eigen::Vector3d &r2) const;

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
  ExclusionList &getExclusions() { return _exclusions; }
  const ExclusionList &getExclusions() const { return _exclusions; }

  BoundaryCondition::eBoxtype getBoxType() const { return _bc->getBoxType(); }

  template <typename iteratable>
  void InsertExclusion(Bead *bead1, iteratable &l);

  bool HasVel() { return _has_vel; }
  void SetHasVel(const bool v) { _has_vel = v; }

  bool HasForce() { return _has_force; }
  void SetHasForce(const bool v) { _has_force = v; }

 protected:
  BoundaryCondition *_bc;

  BoundaryCondition::eBoxtype autoDetectBoxType(
      const Eigen::Matrix3d &box) const;

  /// bead types in the topology
  std::unordered_map<std::string, Index> beadtypes_;

  /// beads in the topology
  BeadContainer _beads;

  /// molecules in the topology
  MoleculeContainer _molecules;

  /// residues in the topology
  ResidueContainer _residues;

  /// bonded interactions in the topology
  InteractionContainer _interactions;

  ExclusionList _exclusions;

  std::map<std::string, Index> _interaction_groups;

  std::map<std::string, std::list<Interaction *> > _interactions_by_group;

  double _time = 0.0;
  Index _step = 0;
  bool _has_vel = false;
  bool _has_force = false;

  /// The particle group (For H5MD file format)
  std::string _particle_group = "unassigned";
};

inline Bead *Topology::CreateBead(Bead::Symmetry symmetry, std::string name,
                                  std::string type, Index resnr, double m,
                                  double q) {

  Bead *b = new Bead(this, _beads.size(), type, symmetry, name, resnr, m, q);
  _beads.push_back(b);
  return b;
}

inline Molecule *Topology::CreateMolecule(std::string name) {
  Molecule *mol = new Molecule(_molecules.size(), name);
  _molecules.push_back(mol);
  return mol;
}

inline Residue *Topology::CreateResidue(std::string name, Index id) {
  Residue *res = new Residue(id, name);
  _residues.push_back(res);
  return res;
}

inline Residue *Topology::CreateResidue(std::string name) {
  Residue *res = new Residue(_residues.size(), name);
  _residues.push_back(res);
  return res;
}

inline Molecule *Topology::MoleculeByIndex(Index index) {
  return _molecules[index];
}

template <typename iteratable>
inline void Topology::InsertExclusion(Bead *bead1, iteratable &l) {
  _exclusions.InsertExclusion(bead1, l);
}

}  // namespace csg
}  // namespace votca

#include "interaction.h"

#endif // VOTCA_CSG_TOPOLOGY_H
