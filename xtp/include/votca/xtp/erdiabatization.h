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
#ifndef VOTCA_XTP_ERDIABATIZATION_H
#define VOTCA_XTP_ERDIABATIZATION_H

#include "logger.h"
#include "votca/xtp/ERIs.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmtool.h"
#include <cstdio>
#include <votca/tools/property.h>
#include <votca/tools/types.h>
#include <votca/xtp/aobasis.h>

namespace votca {
namespace xtp {

class ERDiabatization {
 public:
  ERDiabatization();

  ERDiabatization(Orbitals& orbitals1, Orbitals& orbitals2, Logger* log,
                  Index state_idx_1, Index state_idx_2, std::string qmtype,
                  bool useRI = true)
      : orbitals1_(orbitals1),
        orbitals2_(orbitals2),
        pLog_(log),
        state_idx_1_{state_idx_1},
        state_idx_2_{state_idx_2},
        qmtype_str_(qmtype),
        useRI_(useRI){};

  // Function to set up the ERI matrices
  void setUpMatrices();

  void configure();

  double Calculate_angle() const;
  Eigen::MatrixXd Calculate_diabatic_H(const double angle) const;

 private:
  Orbitals& orbitals1_;
  Orbitals& orbitals2_;
  QMStateType qmtype_;
  Logger* pLog_;

  ERIs eris_;
  AOBasis dftbasis_;
  AOBasis auxbasis_;

  Index state_idx_1_;
  Index state_idx_2_;

  std::string qmtype_str_;

  bool hasRI_;
  bool useRI_;

  double E1_;
  double E2_;

  double CalculateR(const Eigen::MatrixXd& D_JK,
                    const Eigen::MatrixXd& D_LM) const;
  Eigen::MatrixXd CalculateU(const double phi) const;
  Eigen::Tensor<double, 4> CalculateRtensor() const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ERDIABATIZATION_H
