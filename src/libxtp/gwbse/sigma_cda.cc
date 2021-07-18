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

#include <votca/tools/constants.h>
#include <votca/xtp/gw.h>
#include <votca/xtp/sigma_cda.h>

namespace votca {
namespace xtp {

// Prepares the zero and imaginary frequency kappa matrices with
// kappa(omega) = epsilon^-1(omega) - 1 needed in numerical
// integration and for the Gaussian tail
void Sigma_CDA::PrepareScreening() {
  ImaginaryAxisIntegration::options opt;
  opt.homo = opt_.homo;
  opt.order = opt_.order;
  opt.qptotal = qptotal_;
  opt.qpmin = opt_.qpmin;
  opt.rpamax = opt_.rpamax;
  opt.rpamin = opt_.rpamin;
  opt.alpha = opt_.alpha;
  opt.quadrature_scheme = opt_.quadrature_scheme;
  // prepare the zero frequency inverse for Gaussian tail
  kDielMxInv_zero_ =
      rpa_.calculate_epsilon_r(std::complex<double>(0.0, 0.0)).inverse();
  kDielMxInv_zero_.diagonal().array() -= 1.0;
  gq_.configure(opt, rpa_, kDielMxInv_zero_);
}

// This function is used in the calculation of the residues and
// calculates the real part of the dielectric function for a complex
// frequency of the kind omega = delta + i*eta. Instead of explicit
// inversion and multiplication with and Imx vector, a linear system
// is solved.
double Sigma_CDA::CalcDiagContribution(
    const Eigen::MatrixXd::ConstRowXpr& Imx_row, double delta,
    double eta) const {
  std::complex<double> delta_eta(delta, eta);

  Eigen::MatrixXd DielMxInv = rpa_.calculate_epsilon_r(delta_eta);
  Eigen::VectorXd x =
      DielMxInv.partialPivLu().solve(Imx_row.transpose()) - Imx_row.transpose();
  return x.dot(Imx_row.transpose());
}

// Step-function prefactor for the residues
double Sigma_CDA::CalcResiduePrefactor(double e_f, double e_m,
                                       double frequency) const {
  double factor = 0.0;
  double tolerance = 1e-10;
  if (e_f < e_m && e_m < frequency) {
    factor = 1.0;
  } else if (e_f > e_m && e_m > frequency) {
    factor = -1.0;
  } else if (std::abs(e_m - frequency) < tolerance && e_f > e_m) {
    factor = -0.5;
  } else if (std::abs(e_m - frequency) < tolerance && e_f < e_m) {
    factor = 0.5;
  }
  return factor;
}

// Calculates the contribution of residues to the correlation
// part of the self-energy of a fixed gw_level
double Sigma_CDA::CalcResidueContribution(double frequency,
                                          Index gw_level) const {

  const Eigen::VectorXd& rpa_energies = rpa_.getRPAInputEnergies();
  Index rpatotal = rpa_energies.size();
  Index gw_level_offset = gw_level + opt_.qpmin - opt_.rpamin;

  double sigma_c = 0.0;
  double sigma_c_tail = 0.0;
  Index homo = opt_.homo - opt_.rpamin;
  Index lumo = homo + 1;
  double fermi_rpa = (rpa_energies(lumo) + rpa_energies(homo)) / 2.0;
  const Eigen::MatrixXd& Imx = Mmn_[gw_level_offset];

  for (Index i = 0; i < rpatotal; ++i) {
    double delta = rpa_energies(i) - frequency;
    double abs_delta = std::abs(delta);
    double factor = CalcResiduePrefactor(fermi_rpa, rpa_energies(i), frequency);

    // Only considering the terms with a abs(prefactor) > 0.
    // The prefactor can be 1,-1,0.5,-0.5 or 0. We avoid calculating the
    // diagonal contribution if the prefactor is 0. We want to calculate it for
    // all the other cases.
    if (std::abs(factor) > 1e-10) {
      sigma_c +=
          factor * CalcDiagContribution(Imx.row(i), abs_delta, rpa_.getEta());
    }
    // adds the contribution from the Gaussian tail
    if (abs_delta > 1e-10) {
      sigma_c_tail +=
          CalcDiagContributionValue_tail(Imx.row(i), delta, opt_.alpha);
    }
  }
  return sigma_c + sigma_c_tail;
}

// Calculates the correlation part of the self-energy for a fixed
// gw_level and frequence by evaluating the numerical integration
// and residue contributions
double Sigma_CDA::CalcCorrelationDiagElement(Index gw_level,
                                             double frequency) const {

  double sigma_c_residue = CalcResidueContribution(frequency, gw_level);
  double sigma_c_integral = gq_.SigmaGQDiag(frequency, gw_level, rpa_.getEta());
  return sigma_c_residue + sigma_c_integral;
}

// Calculates the contribuion of the tail correction to the
// residue term
double Sigma_CDA::CalcDiagContributionValue_tail(
    const Eigen::MatrixXd::ConstRowXpr& Imx_row, double delta,
    double alpha) const {

  double erfc_factor = 0.5 * std::copysign(1.0, delta) *
                       std::exp(std::pow(alpha * delta, 2)) *
                       std::erfc(std::abs(alpha * delta));

  double value = (Imx_row * kDielMxInv_zero_).dot(Imx_row);
  return value * erfc_factor;
}

}  // namespace xtp
}  // namespace votca
