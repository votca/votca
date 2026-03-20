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
#ifndef VOTCA_XTP_SIGMA_EXACT_UKS_H
#define VOTCA_XTP_SIGMA_EXACT_UKS_H

#include "votca/xtp/rpa_uks.h"
#include "votca/xtp/sigma_base_uks.h"

namespace votca {
namespace xtp {

class Sigma_Exact_UKS : public Sigma_base_UKS {
 public:
  Sigma_Exact_UKS(TCMatrix_gwbse_spin& Mmn, RPA_UKS& rpa,
                  TCMatrix::SpinChannel spin)
      : Sigma_base_UKS(Mmn, rpa, spin) {}

  void PrepareScreening() final;

  double CalcCorrelationDiagElement(Index gw_level,
                                    double frequency) const final;

  double CalcCorrelationDiagElementDerivative(Index gw_level,
                                              double frequency) const final;

  double CalcCorrelationOffDiagElement(Index gw_level1, Index gw_level2,
                                       double frequency1,
                                       double frequency2) const final;

 private:
  // void BuildScreeningModes(const Eigen::MatrixXd& XpY,
  //                          const Eigen::VectorXd& omegas);

  Eigen::VectorXd rpa_omegas_;
  std::vector<Eigen::VectorXd> screening_modes_;
  std::vector<Eigen::MatrixXd> residues_;
};

}  // namespace xtp
}  // namespace votca

#endif