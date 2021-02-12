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
#ifndef VOTCA_CSG_INTERACTION_H
#define VOTCA_CSG_INTERACTION_H

// Standard includes
#include <sstream>
#include <string>

// Local VOTCA includes
#include "bead.h"

namespace votca {
namespace csg {

class Topology;
/**
    \brief base class for all interactions

    This is the base class for all interactions.

    \todo double names/groups right, add molecules!!
*/
class Interaction {
 public:
  virtual double EvaluateVar(const Topology &top) const = 0;

  std::string getName() const { return _name; }

  void setGroup(const std::string &group) {
    _group = group;
    RebuildName();
  }
  const std::string &getGroup() const {
    assert(_group.compare("") != 0);
    return _group;
  }

  // the group id is set by topology, when interaction is added to it
  // \todo if the group name is changed later, group id should be updated by
  // topology
  Index getGroupId() {
    assert(_group_id != -1);
    return _group_id;
  }
  void setGroupId(Index id) { _group_id = id; }

  void setIndex(const Index &index) {
    _index = index;
    RebuildName();
  }
  const Index &getIndex() const {
    assert(_index != -1);
    return _index;
  }

  void setMolecule(const Index &mol) {
    _mol = mol;
    RebuildName();
  }
  const Index &getMolecule() const {
    assert(_mol != -1);
    return _mol;
  }

  virtual Eigen::Vector3d Grad(const Topology &top, Index bead) const = 0;
  Index BeadCount() const { return _beads.size(); }
  Index getBeadId(Index bead) const {
    assert(bead > -1 && boost::lexical_cast<size_t>(bead) < _beads.size());
    return _beads[bead];
  }

 protected:
  template <class T>
  void construct(T &beads) {
    _beads.reserve(beads.size());
    Index i = 0;
    for (const Index &bead_id : beads) {
      _beads[i] = bead_id;
      ++i;
    }
  }
  friend class Topology;

  Index _index = -1;
  std::string _group = "";
  Index _group_id = -1;
  std::string _name = "";
  Index _mol = -1;
  std::vector<Index> _beads;

  void RebuildName();
};

/**
    \brief bond interaction
*/
class IBond : public Interaction {
 public:
  double EvaluateVar(const Topology &top) const override;
  Eigen::Vector3d Grad(const Topology &top, Index bead) const override;

 private:
  template <class T>
  IBond(T &beads) {
    assert(beads.size() == 2);
    construct(beads);
  }
  friend class Topology;
};

/**
    \brief angle interaction
*/
class IAngle : public Interaction {
 public:
  double EvaluateVar(const Topology &top) const override;
  Eigen::Vector3d Grad(const Topology &top, Index bead) const override;

 private:
  template <class T>
  IAngle(T &beads) {
    assert(beads.size() == 3);
    construct(beads);
  }

  friend class Topology;
};

/**
    \brief dihedral interaction
*/
class IDihedral : public Interaction {
 public:
  double EvaluateVar(const Topology &top) const override;
  Eigen::Vector3d Grad(const Topology &top, Index bead) const override;

 private:
  template <class T>
  IDihedral(T &beads) {
    assert(beads.size() == 4);
    construct(beads);
  }

  friend class Topology;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_INTERACTION_H
