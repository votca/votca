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

#include <array>

#include "eigen.h"
#include "votca/xtp/orbitals.h"

namespace votca {
namespace xtp {

typedef std::array<std::vector<std::array<Index, 2>>, 5> OrbTranspositions;

class OrbReorder {
 public:
  OrbReorder(OrbTranspositions transpositions,
             std::array<Index, 25> multipliers)
      : _multipliers(multipliers) {
    _transpositions = transpositions;
  }

  ~OrbReorder() = default;

  void reorderOrbitals(Eigen::MatrixXd& moCoefficients, AOBasis& basis);

 private:
  std::array<Index, 25> _multipliers;
  OrbTranspositions _transpositions;

  std::vector<std::array<Index, 2>> getTranspositions(Index nrOfFunctions);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ORBREODER_H