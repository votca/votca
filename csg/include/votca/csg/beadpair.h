/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_CSG_BEADPAIR_H
#define VOTCA_CSG_BEADPAIR_H

// VOTCA includes
#include <votca/tools/eigen.h>

namespace votca {
namespace csg {

/**
   \brief A particle pair

   This class defines a particle pair. The future plan is, that the Pair class
   can be overloaded and Particle list creates these inherited pairs.

 */

class BeadPair {
 public:
  BeadPair() = default;
  BeadPair(Bead *bead1, Bead *bead2, Eigen::Vector3d r)
      : pair_(std::pair<Bead *, Bead *>(bead1, bead2)),
        r_(r),
        dist_(r.norm()) {}

  Bead *first() { return pair_.first; }
  Bead *second() { return pair_.second; }
  /// \brief the vector connecting two beads
  const Eigen::Vector3d &r() const { return r_; }
  /// \brief the distance of the beads
  double dist() const { return dist_; }

 protected:
  std::pair<Bead *, Bead *> pair_;

  Eigen::Vector3d r_;
  double dist_;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADPAIR_H
