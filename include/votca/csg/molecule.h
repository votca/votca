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

#ifndef VOTCA_CSG_MOLECULE_H
#define VOTCA_CSG_MOLECULE_H

#include "bead.h"
#include <assert.h>
#include <map>
#include <string>
#include <vector>

namespace votca {
namespace csg {

class Interaction;

/**
    \brief information about molecules

    The Molecule class stores which beads belong to a molecule.
    The organization of beads into molecules is needed for the CG mapping.

    \todo sort atoms in molecule

*/
class Molecule : public TopologyItem {
 public:
  /// get the molecule ID
  Index getId() const { return _id; }

  /// get the name of the molecule
  const std::string &getName() const { return _name; }

  /// set the name of the molecule
  void setName(const std::string &name) { _name = name; }

  /// Add a bead to the molecule
  void AddBead(Bead *bead, const std::string &name);
  /// get the id of a bead in the molecule
  Bead *getBead(Index bead) { return _beads[bead]; }
  Index getBeadId(Index bead) { return _beads[bead]->getId(); }
  Index getBeadIdByName(const std::string &name);

  /// get the number of beads in the molecule
  Index BeadCount() const { return _beads.size(); }

  const std::vector<Bead *> &Beads() const { return _beads; }
  std::vector<Bead *> &Beads() { return _beads; }
  /// find a bead by it's name
  Index getBeadByName(const std::string &name);
  std::string getBeadName(Index bead) { return _bead_names[bead]; }

  /// Add an interaction to the molecule
  void AddInteraction(Interaction *ic) { _interactions.push_back(ic); }

  std::vector<Interaction *> Interactions() { return _interactions; }

  template <typename T>
  void setUserData(T *userdata) {
    _userdata = (void *)userdata;
  }

  template <typename T>
  T *getUserData() {
    return (T *)_userdata;
  }

 private:
  // maps a name to a bead id
  std::map<std::string, Index> _beadmap;
  std::vector<Interaction *> _interactions;

  // id of the molecules
  Index _id;

  // name of the molecule
  std::string _name;
  // the beads in the molecule
  std::vector<Bead *> _beads;
  std::vector<std::string> _bead_names;

  void *_userdata;

  /// constructor
  Molecule(Topology *parent, Index id, std::string name)
      : TopologyItem(parent), _id(id), _name(name) {}

  friend class Topology;
};

inline Index Molecule::getBeadIdByName(const std::string &name) {
  Index i = getBeadByName(name);
  if (i < 0) {
    {
      return i;
    }
  }
  return _beads[i]->getId();
}

}  // namespace csg
}  // namespace votca

#endif // VOTCA_CSG_MOLECULE_H 
