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

void Sigma_CDA::PrepareScreening() {
  GaussianQuadrature::options opt;
  opt.homo = _opt.homo;
  opt.order = _opt.order;
  opt.qptotal = _qptotal;
  opt.qpmin = _opt.qpmin;
  opt.rpamin = _opt.rpamin;
  opt.alpha = _opt.alpha;
  opt.quadrature_scheme = _opt.quadrature_scheme;

  // prepare the zero frequency inverse for Gaussian tail
  _kDielMxInv_zero =
      _rpa.calculate_epsilon_r(std::complex<double>(0.0, 0.0)).inverse();
  _kDielMxInv_zero.diagonal().array() -= 1.0;

  _gq.configure(opt, _rpa, _kDielMxInv_zero);
}

// This function is used in the calculation of the residues
// This calculates eps^-1 (inverse of the dielectric function) for complex
// frequencies of the kind omega = delta + i*eta
double Sigma_CDA::CalcDiagContribution(
    const Eigen::MatrixXd::ConstRowXpr& Imx_row, double delta,
    double eta) const {
  std::complex<double> delta_eta(delta, eta);

  Eigen::MatrixXd DielMxInv = _rpa.calculate_epsilon_r(delta_eta);
  // tried linear system solvers
  // - llt() -> fails
  // - ldlt() -> works
  // - ColPivHouseholderQR -> works
  // - HouseholderQR -> works
  // - partialPivLU -> works, should be fastest according to Eigen benchmark
  Eigen::VectorXd x =
      DielMxInv.partialPivLu().solve(Imx_row.transpose()) - Imx_row.transpose();
  return x.dot(Imx_row.transpose());

  // Eigen::MatrixXd DielMxInv = _rpa.calculate_epsilon_r(delta_eta).inverse();
  // DielMxInv.diagonal().array() -= 1.0;
  // return ((Imx_row * DielMxInv).cwiseProduct(Imx_row)).sum();
}

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

double Sigma_CDA::CalcResidueContribution(double frequency,
                                          Index gw_level) const {

  const Eigen::VectorXd& rpa_energies = _rpa.getRPAInputEnergies();
  Index rpatotal = rpa_energies.size();
  Index gw_level_offset = gw_level + _opt.qpmin - _opt.rpamin;

  double sigma_c = 0.0;
  double sigma_c_tail = 0.0;
  Index homo = _opt.homo - _opt.rpamin;
  Index lumo = homo + 1;
  double fermi_rpa = (rpa_energies(lumo) + rpa_energies(homo)) / 2.0;
  const Eigen::MatrixXd& Imx = _Mmn[gw_level_offset];

  for (Index i = 0; i < rpatotal; ++i) {
    double delta = rpa_energies(i) - frequency;
    double abs_delta = std::abs(delta);
    double factor = CalcResiduePrefactor(fermi_rpa, rpa_energies(i), frequency);

    // Only considering the terms with a abs(prefactor) > 0.
    // The prefactor can be 1,-1,0.5,-0.5 or 0. We avoid calculating the
    // diagonal contribution if the prefactor is 0. We want to calculate it for
    // all the other cases.
    if (std::abs(factor) > 1e-10) {
      sigma_c += factor * CalcDiagContribution(Imx.row(i), abs_delta, _eta);
    }
    // This part should allow to add a smooth tail
    if (abs_delta > 1e-10) {  // feeding delta = 0 into copysign gives
                              // interesting results
      sigma_c_tail +=
          CalcDiagContributionValue_tail(Imx.row(i), delta, _opt.alpha);
    }
  }
  return sigma_c + sigma_c_tail;
}

double Sigma_CDA::CalcCorrelationDiagElement(Index gw_level,
                                             double frequency) const {

  double sigma_c_residue = CalcResidueContribution(frequency, gw_level);

  double sigma_c_integral = _gq.SigmaGQDiag(frequency, gw_level, _eta);

  return sigma_c_residue + sigma_c_integral;
}

double Sigma_CDA::CalcDiagContributionValue_tail(Eigen::RowVectorXd Imx_row,
                                                 double delta,
                                                 double alpha) const {

  double erfc_factor = 0.5 * std::copysign(1.0, delta) *
                       std::exp(std::pow(alpha * delta, 2)) *
                       std::erfc(std::abs(alpha * delta));

  double value = ((Imx_row * _kDielMxInv_zero).cwiseProduct(Imx_row)).sum();

  return value * erfc_factor;
}

}  // namespace xtp
}  // namespace votca
