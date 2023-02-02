/*
 *            Copyright 2009-2022 The VOTCA Development Team
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
#ifndef VOTCA_XTP_TRANSITION_DENSITIES_H
#define VOTCA_XTP_TRANSITION_DENSITIES_H

// Standard includes
#include <memory>

// Local VOTCA includes
#include "logger.h"
#include "orbitals.h"
#include "qmstate.h"

namespace votca {
namespace xtp {
/**
 *  \brief  Generalized transition densities tools for different excited states
 *
 *
 */

class TransitionDensities {
 public:
  TransitionDensities();

  TransitionDensities(Orbitals& orbitals1, Orbitals& orbitals2, Logger* log)
      : orbitals1_(orbitals1), orbitals2_(orbitals2), log_(log){};

  void configure();

  Eigen::MatrixXd Matrix(QMState state1, QMState state2);

 private:
  Orbitals& orbitals1_;
  Orbitals& orbitals2_;
  Logger* log_;

  AOBasis dftbasis_;

  Index bse_cmax_;
  Index bse_cmin_;
  Index bse_vmax_;
  Index bse_vmin_;
  Index bse_vtotal_;
  Index bse_ctotal_;
  Index basissize_;

  Eigen::MatrixXd occlevels1_;
  Eigen::MatrixXd virtlevels1_;

  Eigen::MatrixXd occlevels2_;
  Eigen::MatrixXd virtlevels2_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TRANSITION_DENSITIES_H
