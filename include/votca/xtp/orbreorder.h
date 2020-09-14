/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_ORBREODER_H
#define VOTCA_XTP_ORBREODER_H

#include <algorithm>
#include <array>

#include "eigen.h"
#include "votca/xtp/orbitals.h"

namespace votca {
namespace xtp {

// structure to store a transposition
using Transposition = std::pair<Index, Index>;

class OrbReorder {
 public:
  OrbReorder(std::array<Index, 25> reorder, std::array<Index, 25> multipliers,
             bool reverse = false);

  ~OrbReorder() = default;

  void reorderOrbitals(Eigen::MatrixXd& moCoefficients, const AOBasis& basis);
  void reorderRowsAndCols(Eigen::MatrixXd& moCoefficients,
                          const AOBasis& basis);

 private:
  // structure to store the transpositions for the first 5 shell types (i.e.
  // s=0, p, d, f, g=4)
  using OrbTranspositions = std::array<std::vector<Transposition>, 5>;
  std::array<Index, 25> _multipliers;
  std::array<Index, 25> _reorder;
  // clang-format off
  // the ordering of the m quantumnumbers for every shell
  std::array<Index, 25> _votcaOrder={
            0, //s
            0,-1,1, //p
            0,-1,1,-2,2, //d
            0,-1,1,-2,2,-3,3, //f 
            0,-1,1,-2,2,-3,3,-4,4 //g
            };
  // clang-format on
  OrbTranspositions _transpositions;

  std::vector<Transposition> computeTranspositions(
      std::vector<Index> vStart, std::vector<Index> vTarget) const;
  std::vector<Index> copySegment(const std::array<Index, 25>& input,
                                 Index start, Index size) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ORBREODER_H