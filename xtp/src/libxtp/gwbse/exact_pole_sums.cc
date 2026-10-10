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
constexpr Index min_batch_for_split = 24;

// Far poles are summed at Chebyshev nodes of the frequency range [c-h, c+h]
// and interpolated. For a pole at distance d = u h from the centre the
// kernel is analytic inside the Bernstein ellipse with
// rho = u + sqrt(u^2 - 1), and n-node interpolation converges like
// rho^-n. Poles are grouped by u; each group uses enough nodes for
// rho^-n < 1e-17 with a margin of about 1.5 in n:
//   u in [2, 3):  rho >= 3.73, 48 nodes  (rho^-48 ~ 3e-28)
//   u in [3, 5):  rho >= 5.83, 32 nodes  (~ 3e-25)
//   u in [5, 10): rho >= 9.90, 24 nodes  (~ 8e-25)
//   u >= 10:      rho >= 19.9, 16 nodes  (~ 2e-21)
// Poles closer than 2 h are summed at every frequency.
constexpr int n_classes = 4;
constexpr std::array<double, n_classes> class_umin = {2.0, 3.0, 5.0, 10.0};
constexpr std::array<Index, n_classes> class_nodes = {48, 32, 24, 16};
constexpr Index max_nodes = 48;

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

Eigen::VectorXd ValuesDirect(const ExactPoles& poles, const double* residues,
                             const Eigen::VectorXd& w,
                             const double* mode_factors) {
  if (mode_factors == nullptr) {
    return ValuesDirect(poles, residues, w);
  }
  const double eta2 = poles.eta * poles.eta;
  const Index B = poles.nlevels();
  const Index F = w.size();
  Eigen::VectorXd acc = Eigen::VectorXd::Zero(F);
  double* out = acc.data();
  const double* x = w.data();
  for (Index s = 0; s < poles.nmodes(); ++s) {
    const double factor = mode_factors[s];
    if (factor == 0.0) {
      continue;
    }
    const double omega = poles.omegas(s);
    const double* r = residues + B * s;
    for (Index m = 0; m < B; ++m) {
      const double a = factor * r[m] * r[m];
      const double z = poles.energies(m) + PoleSign(m, poles.n_occ) * omega;
#pragma omp simd
      for (Index f = 0; f < F; ++f) {
        out[f] += a * Kernel(x[f] - z, eta2);
      }
    }
  }
  return acc;
}

Eigen::VectorXd Values(const ExactPoles& poles, const double* residues,
                       const Eigen::VectorXd& w) {
  return Values(poles, residues, w, nullptr);
}

Eigen::VectorXd Values(const ExactPoles& poles, const double* residues,
                       const Eigen::VectorXd& w, const double* mode_factors) {
  const Index F = w.size();
  if (F < min_batch_for_split) {
    return ValuesDirect(poles, residues, w, mode_factors);
  }
  const double lo = w.minCoeff();
  const double hi = w.maxCoeff();
  const double c = 0.5 * (lo + hi);
  const double h = 0.5 * (hi - lo);
  if (!(h > 0.0)) {
    return ValuesDirect(poles, residues, w, mode_factors);
  }
  const double eta2 = poles.eta * poles.eta;
  const Index B = poles.nlevels();
  const double inv_h = 1.0 / h;

  // nodes and their angles per class
  std::array<std::array<double, max_nodes>, n_classes> node{};
  std::array<std::array<double, max_nodes>, n_classes> theta{};
  for (int k = 0; k < n_classes; ++k) {
    const Index n = class_nodes[std::size_t(k)];
    for (Index j = 0; j < n; ++j) {
      theta[std::size_t(k)][std::size_t(j)] =
          M_PI * (double(j) + 0.5) / double(n);
      node[std::size_t(k)][std::size_t(j)] =
          c + h * std::cos(theta[std::size_t(k)][std::size_t(j)]);
    }
  }
  std::array<std::array<double, max_nodes>, n_classes> far{};
  Eigen::VectorXd acc = Eigen::VectorXd::Zero(F);
  double* out = acc.data();
  const double* x = w.data();
  for (Index s = 0; s < poles.nmodes(); ++s) {
    const double factor = (mode_factors == nullptr) ? 1.0 : mode_factors[s];
    if (factor == 0.0) {
      continue;
    }
    const double omega = poles.omegas(s);
    const double* r = residues + B * s;
    for (Index m = 0; m < B; ++m) {
      const double a = factor * r[m] * r[m];
      const double z = poles.energies(m) + PoleSign(m, poles.n_occ) * omega;
      const double u = std::abs(z - c) * inv_h;
      if (u < class_umin[0]) {
#pragma omp simd
        for (Index f = 0; f < F; ++f) {
          out[f] += a * Kernel(x[f] - z, eta2);
        }
        continue;
      }
      const int k = (u >= class_umin[3])   ? 3
                    : (u >= class_umin[2]) ? 2
                    : (u >= class_umin[1]) ? 1
                                           : 0;
      const Index n = class_nodes[std::size_t(k)];
      const double* nk = node[std::size_t(k)].data();
      double* fk = far[std::size_t(k)].data();
#pragma omp simd
      for (Index j = 0; j < n; ++j) {
        fk[j] += a * Kernel(nk[j] - z, eta2);
      }
    }
  }
  // per class: Chebyshev coefficients of the far part on [lo, hi], summed
  // into one series of max_nodes terms, evaluated by Clenshaw
  std::array<double, max_nodes> coef{};
  for (int k = 0; k < n_classes; ++k) {
    const Index n = class_nodes[std::size_t(k)];
    for (Index j = 0; j < n; ++j) {
      double sum = 0.0;
      for (Index i = 0; i < n; ++i) {
        sum += far[std::size_t(k)][std::size_t(i)] *
               std::cos(double(j) * theta[std::size_t(k)][std::size_t(i)]);
      }
      coef[std::size_t(j)] += 2.0 * sum / double(n) * ((j == 0) ? 0.5 : 1.0);
    }
  }
  for (Index f = 0; f < F; ++f) {
    const double u = (x[f] - c) * inv_h;
    double b1 = 0.0;
    double b2 = 0.0;
    for (Index j = max_nodes - 1; j >= 1; --j) {
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
