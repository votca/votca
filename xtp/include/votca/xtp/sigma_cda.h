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

#ifndef VOTCA_XTP_SIGMA_CDA_H
#define VOTCA_XTP_SIGMA_CDA_H
#include "votca/xtp/ImaginaryAxisIntegration.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/sigma_base.h"
#include <complex>

// This computes the whole expectation matrix for the correlational part of the
// self-energy with the Contour Deformation Approach according to Eqns 28 and 29
// of JCP 152, 114103 (2020). There are two contributions:
// - from a numerical integration using Gaussian quadratures along the imaginary
//   frequncy axis (Eq. 28)
// - from the residues included in the contours (Eq.29)
// Both contributions contain term from a Gaussian tail with parameter alpha.
namespace votca {
namespace xtp {

class Sigma_CDA : public Sigma_base {

 public:
  Sigma_CDA(TCMatrix_gwbse& Mmn, RPA& rpa)
      : Sigma_base(Mmn, rpa), _gq(rpa.getRPAInputEnergies(), Mmn){};

  ~Sigma_CDA() = default;

  // Prepares the zero and imaginary frequency kappa matrices with
  // kappa(omega) = epsilon^-1(omega) - 1 needed in numerical
  // integration and for the Gaussian tail
  void PrepareScreening() final;

  // calculates the diagonal elements of the self-energy correlation part
  double CalcCorrelationDiagElement(Index gw_level,
                                    double frequency) const final;

  // numerical derivatice of the self-energy
  double CalcCorrelationDiagElementDerivative(Index gw_level,
                                              double frequency) const final {
    double h = 1e-3;
    double plus = CalcCorrelationDiagElement(gw_level, frequency + h);
    double minus = CalcCorrelationDiagElement(gw_level, frequency - h);
    return (plus - minus) / (2 * h);
  }
  // Calculates Sigma_c off-diagonal elements
  double CalcCorrelationOffDiagElement(Index, Index, double,
                                       double) const final {
    return 0;
  }

 private:
  // Theta-function weight of a residue
  double CalcResiduePrefactor(double e_f, double e_m, double frequency) const;

  // Sigma_c from all possible residues for given gw_level and frequency
  double CalcResidueContribution(double frequency, Index gw_level) const;

  // Sigma_c part from a single residue for a given gw_level and frequency
  double CalcDiagContribution(const Eigen::MatrixXd::ConstRowXpr& Imx_row,
                              double delta, double eta) const;

  // Sigma_c part from Gaussian tail correction
  double CalcDiagContributionValue_tail(
      const Eigen::MatrixXd::ConstRowXpr& Imx_row, double delta,
      double alpha) const;

  ImaginaryAxisIntegration _gq;
  Eigen::MatrixXd _kDielMxInv_zero;  // kappa = eps^-1 - 1 matrix
};

}  // namespace xtp

}  // namespace votca

#endif  // VOTCA_XTP_SIGMA_CDA_H