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
#ifndef VOTCA_XTP_GMHDIABATIZATION_H
#define VOTCA_XTP_GMHDIABATIZATION_H

#include "logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmtool.h"
#include <cstdio>
#include <votca/tools/property.h>
#include <votca/tools/types.h>
#include <votca/xtp/aobasis.h>

namespace votca {
namespace xtp {

class GMHDiabatization {
 public:
  GMHDiabatization();

  GMHDiabatization(Orbitals& orbitals1, Orbitals& orbitals2, Logger* log,
                   Index state_idx_1, Index state_idx_2, std::string qmstate)
      : orbitals1_(orbitals1),
        orbitals2_(orbitals2),
        pLog_(log),
        state_idx_1_(state_idx_1),
        state_idx_2_(state_idx_2),
        qmstate_str_(qmstate) {};

  void configure();

  std::pair<double, double> calculate_coupling();
  const std::pair<double, double> adiabatic_energies() {
    return std::pair<double, double>(E1_, E2_);
  };

 private:
  Orbitals& orbitals1_;
  Orbitals& orbitals2_;
  QMStateType qmtype_;
  Logger* pLog_;

  AOBasis dftbasis_;

  Index state_idx_1_;
  Index state_idx_2_;
  std::string qmstate_str_;

  double E1_;
  double E2_;

  Eigen::Vector3d transition_dipole(QMState state1, QMState state2);
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GMHDIABATIZATION_H
