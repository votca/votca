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

// Standard includes
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

// Local VOTCA includes
#include "votca/xtp/exact_pole_sums.h"

namespace votca {
namespace xtp {
namespace exact_pole_sums {

namespace {

// Batches smaller than this are summed directly
constexpr Index min_batch_for_split = 64;
// Chebyshev nodes for the far poles. Their singularities (z +- i eta) lie at
// least one half-width outside the range, i.e. outside the Bernstein ellipse
// with rho = 2 + sqrt(3); with rho = 3 the error is below 3^-48 ~ 1e-23 of
// the far sum's magnitude on that ellipse.
constexpr Index cheb_nodes = 48;

inline double Kernel(double t, double eta2) { return t / (t * t + eta2); }

// sign of Omega in the pole position: -1 occupied, +1 virtual
inline double PoleSign(Index m, Index n_occ) {
  return (m < n_occ) ? -1.0 : 1.0;
}

}  // namespace

double Value(const ExactPoles& poles, const double* residues, double w) {
  const double eta2 = poles.eta * poles.eta;
  const Index B = poles.nlevels();
  // t = (w - e_m) - sign_m Omega_s
  std::vector<double> d(static_cast<std::size_t>(B));
  std::vector<double> sgn(static_cast<std::size_t>(B));
  for (Index m = 0; m < B; ++m) {
    d[std::size_t(m)] = w - poles.energies(m);
    sgn[std::size_t(m)] = PoleSign(m, poles.n_occ);
  }
  double sum = 0.0;
  for (Index s = 0; s < poles.nmodes(); ++s) {
    const double omega = poles.omegas(s);
    const double* r = residues + B * s;
    double part = 0.0;
#pragma omp simd reduction(+ : part)
    for (Index m = 0; m < B; ++m) {
      const double t = d[std::size_t(m)] - sgn[std::size_t(m)] * omega;
      part += r[m] * r[m] * Kernel(t, eta2);
    }
    sum += part;
  }
  return sum;
}

double Derivative(const ExactPoles& poles, const double* residues, double w) {
  const double eta2 = poles.eta * poles.eta;
  const Index B = poles.nlevels();
  std::vector<double> d(static_cast<std::size_t>(B));
  std::vector<double> sgn(static_cast<std::size_t>(B));
  for (Index m = 0; m < B; ++m) {
    d[std::size_t(m)] = w - poles.energies(m);
    sgn[std::size_t(m)] = PoleSign(m, poles.n_occ);
  }
  double sum = 0.0;
  for (Index s = 0; s < poles.nmodes(); ++s) {
    const double omega = poles.omegas(s);
    const double* r = residues + B * s;
    double part = 0.0;
#pragma omp simd reduction(+ : part)
    for (Index m = 0; m < B; ++m) {
      const double t = d[std::size_t(m)] - sgn[std::size_t(m)] * omega;
      const double t2 = t * t;
      const double denom = t2 + eta2;
      part += r[m] * r[m] * (eta2 - t2) / (denom * denom);
    }
    sum += part;
  }
  return sum;
}

namespace {
// sum_p weight(r_p) f(w - z_p) at every frequency
template <class Weight>
Eigen::VectorXd DirectSum(const ExactPoles& poles, const double* r_all,
                          const Eigen::VectorXd& w, Weight weight) {
  const double eta2 = poles.eta * poles.eta;
  const Index B = poles.nlevels();
  const Index F = w.size();
  Eigen::VectorXd acc = Eigen::VectorXd::Zero(F);
  double* out = acc.data();
  const double* x = w.data();
  for (Index s = 0; s < poles.nmodes(); ++s) {
    const double omega = poles.omegas(s);
    const double* r = r_all + B * s;
    for (Index m = 0; m < B; ++m) {
      const double a = weight(r[m]);
      const double z = poles.energies(m) + PoleSign(m, poles.n_occ) * omega;
#pragma omp simd
      for (Index f = 0; f < F; ++f) {
        out[f] += a * Kernel(x[f] - z, eta2);
      }
    }
  }
  return acc;
}
}  // namespace

Eigen::VectorXd ValuesDirect(const ExactPoles& poles, const double* residues,
                             const Eigen::VectorXd& w) {
  return DirectSum(poles, residues, w, [](double r) { return r * r; });
}

Eigen::VectorXd ValuesDirectWeighted(const ExactPoles& poles,
                                     const double* weights,
                                     const Eigen::VectorXd& w) {
  return DirectSum(poles, weights, w, [](double a) { return a; });
}

Eigen::VectorXd Values(const ExactPoles& poles, const double* residues,
                       const Eigen::VectorXd& w) {
  const Index F = w.size();
  if (F < min_batch_for_split) {
    return ValuesDirect(poles, residues, w);
  }
  const double lo = w.minCoeff();
  const double hi = w.maxCoeff();
  const double c = 0.5 * (lo + hi);
  const double h = 0.5 * (hi - lo);
  if (!(h > 0.0)) {
    return ValuesDirect(poles, residues, w);
  }
  const double near = 2.0 * h;  // |z - c| below this: near pole
  const double eta2 = poles.eta * poles.eta;
  const Index B = poles.nlevels();

  std::array<double, cheb_nodes> node;
  std::array<double, cheb_nodes> theta;
  for (Index k = 0; k < cheb_nodes; ++k) {
    theta[std::size_t(k)] = M_PI * (double(k) + 0.5) / double(cheb_nodes);
    node[std::size_t(k)] = c + h * std::cos(theta[std::size_t(k)]);
  }
  std::array<double, cheb_nodes> far{};
  Eigen::VectorXd acc = Eigen::VectorXd::Zero(F);
  double* out = acc.data();
  const double* x = w.data();
  for (Index s = 0; s < poles.nmodes(); ++s) {
    const double omega = poles.omegas(s);
    const double* r = residues + B * s;
    for (Index m = 0; m < B; ++m) {
      const double a = r[m] * r[m];
      const double z = poles.energies(m) + PoleSign(m, poles.n_occ) * omega;
      if (std::abs(z - c) < near) {
#pragma omp simd
        for (Index f = 0; f < F; ++f) {
          out[f] += a * Kernel(x[f] - z, eta2);
        }
      } else {
#pragma omp simd
        for (Index k = 0; k < cheb_nodes; ++k) {
          far[std::size_t(k)] += a * Kernel(node[std::size_t(k)] - z, eta2);
        }
      }
    }
  }
  // Chebyshev coefficients of the far part on [lo, hi]
  std::array<double, cheb_nodes> coef{};
  for (Index j = 0; j < cheb_nodes; ++j) {
    double sum = 0.0;
    for (Index k = 0; k < cheb_nodes; ++k) {
      sum += far[std::size_t(k)] * std::cos(double(j) * theta[std::size_t(k)]);
    }
    coef[std::size_t(j)] = 2.0 * sum / double(cheb_nodes);
  }
  coef[0] *= 0.5;
  // Clenshaw
  for (Index f = 0; f < F; ++f) {
    const double u = (x[f] - c) / h;
    double b1 = 0.0;
    double b2 = 0.0;
    for (Index j = cheb_nodes - 1; j >= 1; --j) {
      const double b0 = 2.0 * u * b1 - b2 + coef[std::size_t(j)];
      b2 = b1;
      b1 = b0;
    }
    out[f] += u * b1 - b2 + coef[0];
  }
  return acc;
}

Eigen::MatrixXd OffDiagonal(const ExactPoles& poles,
                            const Eigen::MatrixXd& residues_all,
                            const Eigen::VectorXd& w) {
  const Index P = residues_all.rows();
  const Index Q = residues_all.cols();
  const Index B = poles.nlevels();
  const double eta2 = poles.eta * poles.eta;
  // chunks of poles of about 32 MB of weighted residues each
  const Index chunk =
      std::max<Index>(256, Index(4e6 / double(std::max<Index>(Q, 1))));
  const Index nchunks = (P + chunk - 1) / chunk;
  Eigen::MatrixXd T = Eigen::MatrixXd::Zero(Q, Q);
#pragma omp parallel
  {
    Eigen::MatrixXd T_local = Eigen::MatrixXd::Zero(Q, Q);
    Eigen::MatrixXd weighted;
#pragma omp for schedule(dynamic)
    for (Index ic = 0; ic < nchunks; ++ic) {
      const Index p0 = ic * chunk;
      const Index np = std::min(chunk, P - p0);
      // weighted(p, b) = r_b(p) f(w_b - z_p)
      weighted.resize(np, Q);
      for (Index b = 0; b < Q; ++b) {
        const double wb = w(b);
        for (Index i = 0; i < np; ++i) {
          const Index p = p0 + i;
          const Index m = p % B;
          const Index s = p / B;
          const double z =
              poles.energies(m) + PoleSign(m, poles.n_occ) * poles.omegas(s);
          weighted(i, b) = residues_all(p, b) * Kernel(wb - z, eta2);
        }
      }
      // T(a, b) += sum_p r_a(p) r_b(p) f(w_b - z_p)
      T_local.noalias() +=
          residues_all.middleRows(p0, np).transpose() * weighted;
    }
#pragma omp critical
    T += T_local;
  }
  Eigen::MatrixXd result = T + T.transpose();
  result.diagonal().setZero();
  return result;
}

}  // namespace exact_pole_sums
}  // namespace xtp
}  // namespace votca
