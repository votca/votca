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

#ifndef VOTCA_XTP_SIGMA_EXACT_H
#define VOTCA_XTP_SIGMA_EXACT_H

// Local VOTCA includes
#include "votca/xtp/exact_pole_sums.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/sigma_base.h"

namespace votca {
namespace xtp {

class TCMatrix_gwbse;
class RPA;

class Sigma_Exact : public Sigma_base {

 public:
  Sigma_Exact(TCMatrix_gwbse& Mmn, RPA& rpa) : Sigma_base(Mmn, rpa){};

  // Sets up the screening parametrisation
  void PrepareScreening() final;
  // Calculates Sigma_c diagonal elements
  double CalcCorrelationDiagElement(Index gw_level,
                                    double frequency) const final;

  double CalcCorrelationDiagElementDerivative(Index gw_level,
                                              double frequency) const final;
  /// All frequencies of a level in one pass (far poles via Chebyshev
  /// interpolation, see exact_pole_sums)
  Eigen::VectorXd CalcCorrelationDiagElements(
      Index gw_level, const Eigen::VectorXd& frequencies) const final;
  // Calculates Sigma_c off-diagonal elements
  double CalcCorrelationOffDiagElement(Index gw_level1, Index gw_level2,
                                       double frequency1,
                                       double frequency2) const final;
  /// All off-diagonal elements as blocked matrix products
  Eigen::MatrixXd CalcCorrelationOffDiag(
      const Eigen::VectorXd& frequencies) const final;

 private:
  Eigen::VectorXd rpa_omegas_;  // Eigenvalues from RPA
  // residues of GW level n (column n): the RPA levels x modes array
  // M_n Z, column-major (pole m + B s)
  Eigen::MatrixXd residues_;

  ExactPoles Poles() const;
  const double* Residues(Index gw_level) const {
    return residues_.col(gw_level).data();
  }

  // Z = sum_v M_v,virt^T (X+Y)_v (aux x RPA modes); the residues of level
  // m are M_m Z
  Eigen::MatrixXd ScreeningModes(const Eigen::MatrixXd& XpY) const;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SIGMA_EXACT_H
