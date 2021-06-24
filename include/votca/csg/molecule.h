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

#pragma once
#ifndef VOTCA_CSG_MOLECULE_H
#define VOTCA_CSG_MOLECULE_H

// Standard includes
#include <cassert>
#include <map>
#include <string>
#include <vector>

// Local VOTCA includes
#include "bead.h"

namespace votca {
namespace csg {

class Interaction;

/**
    \brief information about molecules

    The Molecule class stores which beads belong to a molecule.
    The organization of beads into molecules is needed for the CG mapping.

    \todo sort atoms in molecule

*/
class Molecule {
 public:
  /// get the molecule ID
  Index getId() const { return id_; }

  /// get the name of the molecule
  const std::string &getName() const { return name_; }

  /// set the name of the molecule
  void setName(const std::string &name) { name_ = name; }

  /// Add a bead to the molecule
  void AddBead(Bead *bead, const std::string &name);
  /// get the id of a bead in the molecule
  Bead *getBead(Index bead) { return beads_[bead]; }
  const Bead *getBead(Index bead) const { return beads_[bead]; }
  Index getBeadId(Index bead) const { return beads_[bead]->getId(); }
  Index getBeadIdByName(const std::string &name) const;

  /// get the number of beads in the molecule
  Index BeadCount() const { return beads_.size(); }

  const std::vector<Bead *> &Beads() const { return beads_; }
  std::vector<Bead *> &Beads() { return beads_; }
  /// find a bead by it's name
  Index getBeadByName(const std::string &name) const;
  std::string getBeadName(const Index bead) { return bead_names_[bead]; }
  const std::string &getBeadName(const Index bead) const {
    return bead_names_[bead];
  }

  /// Add an interaction to the molecule
  void AddInteraction(Interaction *ic) { interactions_.push_back(ic); }

  std::vector<Interaction *> Interactions() { return interactions_; }
  const std::vector<Interaction *> &Interactions() const {
    return interactions_;
  }

 private:
  // maps a name to a bead id
  std::map<std::string, Index> beadmap_;
  std::vector<Interaction *> interactions_;

  // id of the molecules
  Index id_;

  // name of the molecule
  std::string name_;
  // the beads in the molecule
  std::vector<Bead *> beads_;
  std::vector<std::string> bead_names_;

  /// constructor
  Molecule(Index id, std::string name) : id_(id), name_(name) {}

  friend class Topology;
};

inline Index Molecule::getBeadIdByName(const std::string &name) const {
  Index i = getBeadByName(name);
  if (i < 0) {
    return i;
  }
  return beads_[i]->getId();
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_MOLECULE_H
