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
#ifndef VOTCA_CSG_BEADLIST_H
#define VOTCA_CSG_BEADLIST_H

// Standard includes
#include <string>
#include <vector>

// VOTCA includes
#include <votca/tools/eigen.h>
#include <votca/tools/types.h>

namespace votca {
namespace csg {
/**
    \brief Generate lists of beads

    This class generates a list of beads based on some criteria, currently
    only the bead type.

*/

class Topology;
class Bead;

class BeadList {
 public:
  /// \brief Select all beads of type "select"
  Index Generate(Topology &top, const std::string &select);
  /// \brief Select all beads of type "select" withn a radius "radius" of
  /// reference vector "ref"
  Index GenerateInSphericalSubvolume(Topology &top, const std::string &select,
                                     Eigen::Vector3d ref, double radius);

  Index size() const { return beads_.size(); }

  bool empty() const { return beads_.empty(); }

  void push_back(Bead *bead) { beads_.push_back(bead); }

  using iterator = typename std::vector<Bead *>::iterator;

  iterator begin() { return beads_.begin(); }
  iterator end() { return beads_.end(); }

  const Topology &getTopology() const { return *topology_; }

 private:
  std::vector<Bead *> beads_;
  Topology *topology_ = nullptr;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADLIST_H
