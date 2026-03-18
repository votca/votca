/*
 *            Copyright 2009-2026 The VOTCA Development Team
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

#include "sigma_cda_uks.h"
#include "votca/xtp/gw.h"

#include <algorithm>
#include <iostream>
#include <limits>

namespace votca {
namespace xtp {

namespace {

constexpr double kResidueFactorTol = 1e-10;
constexpr double kCdaWarnMinDelta = 1e-4;  // Ha
constexpr double kCdaWarnKappa = 1e3;
constexpr double kTiny = 1e-14;

const char* SpinName(TCMatrix::SpinChannel spin) {
  return (spin == TCMatrix::SpinChannel::Alpha) ? "alpha" : "beta";
}

}  // namespace

void Sigma_CDA_UKS::PrepareScreening() {
  ImaginaryAxisIntegration::options opt;
  opt.homo = opt_.homo;
  opt.order = opt_.order;
  opt.qptotal = qptotal_;
  opt.qpmin = opt_.qpmin;
  opt.rpamax = opt_.rpamax;
  opt.rpamin = opt_.rpamin;
  opt.alpha = opt_.alpha;
  opt.quadrature_scheme = opt_.quadrature_scheme;

  kDielMxInv_zero_ =
      rpa_.calculate_epsilon_r(std::complex<double>(0.0, 0.0)).inverse();
  kDielMxInv_zero_.diagonal().array() -= 1.0;
  gq_.configure(opt, rpa_, kDielMxInv_zero_);
}

double Sigma_CDA_UKS::CalcDiagContribution(
    const Eigen::MatrixXd::ConstRowXpr& Imx_row, double delta,
    double eta) const {
  std::complex<double> delta_eta(delta, eta);

  Eigen::MatrixXd DielMxInv = rpa_.calculate_epsilon_r(delta_eta);
  Eigen::VectorXd x =
      DielMxInv.partialPivLu().solve(Imx_row.transpose()) - Imx_row.transpose();

  return x.dot(Imx_row.transpose());
}

double Sigma_CDA_UKS::CalcResiduePrefactor(double e_f, double e_m,
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

double Sigma_CDA_UKS::CalcResidueContribution(double frequency,
                                              Index gw_level) const {
  const Eigen::VectorXd& rpa_energies = getSpinRPAInputEnergies();
  Index rpatotal = rpa_energies.size();
  Index gw_level_offset = gw_level + opt_.qpmin - opt_.rpamin;

  double sigma_c = 0.0;
  double sigma_c_tail = 0.0;
  Index homo = opt_.homo - opt_.rpamin;
  Index lumo = homo + 1;
  double fermi_rpa = (rpa_energies(lumo) + rpa_energies(homo)) / 2.0;

  double min_abs_delta = std::numeric_limits<double>::infinity();
  double sum_abs_contributions = 0.0;
  double sum_signed_contributions = 0.0;

  const Eigen::MatrixXd& Imx = Mmn_[gw_level_offset];

  for (Index i = 0; i < rpatotal; ++i) {
    double delta = rpa_energies(i) - frequency;
    double abs_delta = std::abs(delta);
    min_abs_delta = std::min(min_abs_delta, abs_delta);

    double factor = CalcResiduePrefactor(fermi_rpa, rpa_energies(i), frequency);

    if (std::abs(factor) > kResidueFactorTol) {
      double diag =
          CalcDiagContribution(Imx.row(i), abs_delta, rpa_.getEta());
      double contribution = factor * diag;
      sigma_c += contribution;

      sum_abs_contributions += std::abs(contribution);
      sum_signed_contributions += contribution;
    }

    if (abs_delta > kResidueFactorTol) {
      double tail =
          CalcDiagContributionValue_tail(Imx.row(i), delta, opt_.alpha);
      sigma_c_tail += tail;
    }
  }

  double kappa = sum_abs_contributions /
                 (std::abs(sum_signed_contributions) + kTiny);

  if (min_abs_delta < kCdaWarnMinDelta || kappa > kCdaWarnKappa) {
    std::cout << "\nWarning: CDA may be unreliable for "
              << SpinName(spin_) << " GW level " << (gw_level + opt_.qpmin)
              << " at omega = " << frequency
              << " Ha: min |e_m - omega| = " << min_abs_delta
              << " Ha, cancellation metric = " << kappa << std::flush;
  }

  return sigma_c + sigma_c_tail;
}

double Sigma_CDA_UKS::CalcCorrelationDiagElement(Index gw_level,
                                                 double frequency) const {
  double sigma_c_residue = CalcResidueContribution(frequency, gw_level);
  double sigma_c_integral = gq_.SigmaGQDiag(frequency, gw_level, rpa_.getEta());
  return sigma_c_residue + sigma_c_integral;
}

double Sigma_CDA_UKS::CalcDiagContributionValue_tail(
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