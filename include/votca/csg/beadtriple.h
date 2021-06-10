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

#ifndef VOTCA_CSG_BEADTRIPLE_H
#define VOTCA_CSG_BEADTRIPLE_H

// Standard includes
#include <tuple>

namespace votca {
namespace csg {

/**
   \brief A triplet of tree Beads

 */

class BeadTriple : public std::tuple<Bead *, Bead *, Bead *> {
 public:
  BeadTriple() = default;
  BeadTriple(Bead *bead1, Bead *bead2, Bead *bead3, Eigen::Vector3d r12,
             Eigen::Vector3d r13, Eigen::Vector3d r23)
      : std::tuple<Bead *, Bead *, Bead *>(bead1, bead2, bead3),
        r12_(r12),
        r13_(r13),
        r23_(r23),
        dist12_(r12.norm()),
        dist13_(r13.norm()),
        dist23_(r23.norm()) {}

  virtual ~BeadTriple() = default;

  /// \brief return the beads
  const Bead *bead1() { return std::get<0>(*this); }
  const Bead *bead2() { return std::get<1>(*this); }
  const Bead *bead3() { return std::get<2>(*this); }

  /// \brief the vector connecting two beads
  Eigen::Vector3d &r12() { return r12_; }
  Eigen::Vector3d &r13() { return r13_; }
  Eigen::Vector3d &r23() { return r23_; }
  /// \brief the distance of the beads
  double &dist12() { return dist12_; }
  double &dist13() { return dist13_; }
  double &dist23() { return dist23_; }

 protected:
  Eigen::Vector3d r12_;
  Eigen::Vector3d r13_;
  Eigen::Vector3d r23_;
  double dist12_;
  double dist13_;
  double dist23_;
};

}  // namespace csg
}  // namespace votca

#endif /*  VOTCA_CSG_BEADTRIPLE_H */
