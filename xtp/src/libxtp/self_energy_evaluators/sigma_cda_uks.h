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
#ifndef VOTCA_XTP_SIGMA_CDA_UKS_H
#define VOTCA_XTP_SIGMA_CDA_UKS_H

#include "votca/xtp/ImaginaryAxisIntegration.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/rpa_uks.h"
#include "votca/xtp/sigma_base_uks.h"
#include <complex>

namespace votca {
namespace xtp {

class Sigma_CDA_UKS : public Sigma_base_UKS {
 public:
  Sigma_CDA_UKS(TCMatrix_gwbse_spin& Mmn, RPA_UKS& rpa,
                TCMatrix::SpinChannel spin)
      : Sigma_base_UKS(Mmn, rpa, spin),
        gq_(getSpinRPAInputEnergies(), Mmn[spin]) {}

  ~Sigma_CDA_UKS() override = default;

  void PrepareScreening() final;

  double CalcCorrelationDiagElement(Index gw_level,
                                    double frequency) const final;

  double CalcCorrelationDiagElementDerivative(Index gw_level,
                                              double frequency) const final {
    double h = 1e-3;
    double plus = CalcCorrelationDiagElement(gw_level, frequency + h);
    double minus = CalcCorrelationDiagElement(gw_level, frequency - h);
    return (plus - minus) / (2 * h);
  }

  double CalcCorrelationOffDiagElement(Index, Index, double,
                                       double) const final {
    return 0.0;
  }

 private:
  double CalcResiduePrefactor(double e_f, double e_m, double frequency) const;
  double CalcResidueContribution(double frequency, Index gw_level) const;
  double CalcDiagContribution(const Eigen::MatrixXd::ConstRowXpr& Imx_row,
                              double delta, double eta,
                              Index gw_level) const;
  double CalcDiagContributionValue_tail(
      const Eigen::MatrixXd::ConstRowXpr& Imx_row, double delta,
      double alpha) const;

  ImaginaryAxisIntegration gq_;
  Eigen::MatrixXd kDielMxInv_zero_;
};

}  // namespace xtp
}  // namespace votca

#endif