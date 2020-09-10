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

#ifndef _VOTCA_XTP_SIGMA_CDA_H
#define _VOTCA_XTP_SIGMA_CDA_H
#include "votca/xtp/gaussian_quadrature.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/sigma_base.h"
#include <complex>

// This computes the whole expectation matrix for the correlational part of the
// self-energy: so, both the residual and the Gauss-Hermite quadrature
// contribution
namespace votca {
namespace xtp {

class TCMatrix_gwbse;
class RPA;

class Sigma_CDA : public Sigma_base {

 public:
  Sigma_CDA(TCMatrix_gwbse& Mmn, RPA& rpa)
      : Sigma_base(Mmn, rpa),
        _gq(rpa.getRPAInputEnergies(), Mmn),
        _eta(rpa.getEta()){};

  ~Sigma_CDA(){};

  void PrepareScreening() final;

  double CalcCorrelationDiagElement(Index gw_level,
                                    double frequency) const final;

  double CalcCorrelationDiagElementDerivative(Index gw_level,
                                              double frequency) const {
    double h = 1e-3;
    double plus = CalcCorrelationDiagElement(gw_level, frequency + h);
    double minus = CalcCorrelationDiagElement(gw_level, frequency - h);
    return (plus - minus) / (2 * h);
  }
  // Calculates Sigma_c off-diagonal elements
  double CalcCorrelationOffDiagElement(Index, Index, double, double) const {
    return 0;
  }

 private:
  double CalcResiduePrefactor(double e_f, double e_m, double frequency) const;

  double CalcResidueContribution( double frequency,
                                 Index gw_level) const;

  double CalcDiagContribution(const Eigen::Ref<const Eigen::MatrixXd>& Imx_row, double delta,
                              double eta) const;

  double CalcDiagContributionValue_tail(Eigen::RowVectorXd Imx_row,
                                                 double delta,
                                                 double alpha) const;
  GaussianQuadrature _gq;

  double _eta;

  Eigen::MatrixXd _kDielMxInv_zero; // kappa = eps^-1 - 1 matrix

};

}  // namespace xtp

}  // namespace votca

#endif /* _VOTCA_XTP_SIGMA_CDA_H */