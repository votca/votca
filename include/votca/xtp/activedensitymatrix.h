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
#ifndef VOTCA_XTP_ACTIVEDENSITYMATRIX_H
#define VOTCA_XTP_ACTIVEDENSITYMATRIX_H
#include "logger.h"
#include "votca/xtp/orbitals.h"

namespace votca {
namespace xtp {
class ActiveDensityMatrix {
 public:
  ActiveDensityMatrix(Orbitals &orbitals, std::vector<Index> activeatoms,
                      Logger &log)
      : orbitals(orbitals), log(log), activeatoms(activeatoms){};
  Eigen::MatrixXd activedensitymatrix(Eigen::MatrixXd &new_mo_coeff);
 Eigen::MatrixXd compute_Dmat_A();

 private:
  Orbitals orbitals;
  Logger &log;
  std::vector<Index> activeatoms;
  BasisSet basis;
  AOBasis aobasis;
};
}  // namespace xtp
}  // namespace votca
#endif