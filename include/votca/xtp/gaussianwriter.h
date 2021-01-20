/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
#ifndef VOTCA_XTP_GAUSSIANWRITER_H
#define VOTCA_XTP_GAUSSIANWRITER_H

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmtool.h"
#include <fstream>

namespace votca {
namespace xtp {
class GaussianWriter {
 public:
  GaussianWriter(Logger& log) : _log(log) { gaussianMultipliers.fill(-1); }

  ~GaussianWriter() = default;

  void WriteFile(const std::string& basename, const Orbitals& orbitals,
                 const QMState state = QMState(QMStateType::statetype::Gstate,
                                               0, false)) const;

 private:
  Logger& _log;
  Index toGaussianL(L l) const;
  std::string reorderedMOCoefficients(const Orbitals& orbitals) const;
  std::string densityMatrixToString(const Orbitals& orbitals,
                                    const QMState& state) const;
  // Setup the reordering parameters
  std::array<Index, 49> gaussianMultipliers;
  // clang-format off
  // the ordering of the m quantumnumbers for every shell in gaussian
  std::array<Index, 49> gaussianOrder={{
            0, //s
            1,-1,0, //p
            0,1,-1,2,-2, //d
            0,1,-1,2,-2,3,-3, //f 
            0,1,-1,2,-2,3,-3,4,-4, //g
            0,1,-1,2,-2,3,-3,4,-4,5,-5, // h
            0,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6 // i
  }};
  // clang-format on
};

}  // namespace xtp
}  // namespace votca

#endif