#include "sigma_ri_reduced.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/threecenter.h"
#include <cmath>
#include <stdexcept>

namespace votca {
namespace xtp {

Sigma_RI_Reduced::Sigma_RI_Reduced(TCMatrix_gwbse& Mmn, RPA& rpa)
    : Sigma_base(Mmn, rpa), rpa_red_(Mmn) {}

std::vector<double> Sigma_RI_Reduced::ImagFrequencyGrid() const {
  std::vector<double> grid;

  if (opt_red_.imag_omega_points <= 0) {
    return grid;
  }

  const double wmax = opt_red_.imag_omega_max;
  const Index N = opt_red_.imag_omega_points;

  grid.reserve(static_cast<size_t>(N));

  for (Index i = 0; i < N; ++i) {
    double x = static_cast<double>(i) / static_cast<double>(N - 1);
    grid.push_back(wmax * x);
  }

  return grid;
}

double Sigma_RI_Reduced::CalcCorrelationDiagElementDerivative(
    Index gw_level, double frequency) const {
  // symmetric finite-difference derivative dSigma/dw
  // step in Hartree; small, but not too small for numerical stability
  const double h = 1e-4;

  const double sigma_p = CalcCorrelationDiagElement(gw_level, frequency + h);
  const double sigma_m = CalcCorrelationDiagElement(gw_level, frequency - h);

  return (sigma_p - sigma_m) / (2.0 * h);
}

std::vector<double> Sigma_RI_Reduced::ImagFrequencyWeights() const {
  const Index N = opt_red_.imag_omega_points;
  std::vector<double> w(N, 1.0);

  if (N <= 1) return w;

  // simple trapezoidal rule
  w.front() *= 0.5;
  w.back() *= 0.5;

  double step = opt_red_.imag_omega_max / static_cast<double>(N - 1);
  for (auto& wi : w) wi *= step;

  return w;
}

void Sigma_RI_Reduced::PrepareScreening() {
  // --- configure reduced RPA ---
  rpa_red_.configure(opt_.homo, opt_.rpamin, opt_.rpamax);
  rpa_red_.setRPAInputEnergies(rpa_.getRPAInputEnergies());

  RPA_RI_Reduced::options ropt;
  ropt.imag_omega_max = opt_red_.imag_omega_max;
  ropt.imag_omega_points = opt_red_.imag_omega_points;
  ropt.basis_threshold = opt_red_.basis_threshold;
  ropt.max_rank = opt_red_.max_rank;

  rpa_red_.configure_reduced(ropt);

  // --- build reduced basis ---
  rpa_red_.BuildReducedBasis();

  // --- build couplings C_im^α ---
  const Index qpoffset = opt_.qpmin - opt_.rpamin;

  reduced_couplings_.resize(qptotal_);

#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; ++gw_level) {
    reduced_couplings_[gw_level] =
        Mmn_[gw_level + qpoffset] * rpa_red_.BasisU();
  }
}

double Sigma_RI_Reduced::CalcCorrelationDiagElement(Index gw_level,
                                                    double frequency) const {
  const Eigen::MatrixXd& Cim = reduced_couplings_[gw_level];
  const Eigen::VectorXd& eps = rpa_red_.getRPAInputEnergies();

  const std::vector<double> omegas = ImagFrequencyGrid();
  const std::vector<double> weights = ImagFrequencyWeights();

  double sigma = 0.0;

  for (Index iw = 0; iw < static_cast<Index>(omegas.size()); ++iw) {
    const double omega = omegas[iw];
    const double w = weights[iw];

    const Eigen::MatrixXd Wc = rpa_red_.BuildReducedWcImag(omega);

    for (Index m = 0; m < rpatotal_; ++m) {
      const Eigen::RowVectorXd c = Cim.row(m);

      const double screened = (c * Wc * c.transpose())(0, 0);

      const double denom = frequency - eps(m);

      // simple imaginary-axis kernel approximation
      sigma += w * screened * denom / (denom * denom + omega * omega);
    }
  }

  return -sigma / M_PI;
}

}  // namespace xtp
}  // namespace votca