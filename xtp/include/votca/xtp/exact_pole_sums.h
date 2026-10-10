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

#pragma once
#ifndef VOTCA_XTP_EXACT_POLE_SUMS_H
#define VOTCA_XTP_EXACT_POLE_SUMS_H

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief Pole sums of the exact (RPA spectral) correlation self-energy.
 *
 * For a GW level n with residues r_n(m, s) (RPA level m, screening mode s),
 *
 *   S_n(w) = sum_{m,s} r_n(m,s)^2 t / (t^2 + eta^2),   t = w - z_{ms},
 *
 * with the pole z_{ms} = e_m - Omega_s for occupied m (m < n_occ) and
 * e_m + Omega_s for virtual m. The residues of one level are a B x S
 * column-major array (B RPA levels, S modes), pole p = m + B s.
 */
struct ExactPoles {
  const Eigen::VectorXd& energies;  // e_m, size B
  Index n_occ;
  const Eigen::VectorXd& omegas;  // Omega_s, size S
  double eta;

  Index nlevels() const { return energies.size(); }
  Index nmodes() const { return omegas.size(); }
  Index npoles() const { return nlevels() * nmodes(); }
};

namespace exact_pole_sums {

/// S_n(w) at one frequency
double Value(const ExactPoles& poles, const double* residues, double w);

/// dS_n/dw at one frequency
double Derivative(const ExactPoles& poles, const double* residues, double w);

/// S_n at many frequencies. For larger batches the poles far from the
/// frequencies (more than the half-width of their range away from it) are
/// summed at Chebyshev nodes of the range and interpolated, which converges
/// exponentially (error far below 1e-14 relative); the near poles are summed
/// at every frequency.
Eigen::VectorXd Values(const ExactPoles& poles, const double* residues,
                       const Eigen::VectorXd& w);

/// The same sum term by term at every frequency (reference for tests)
Eigen::VectorXd ValuesDirect(const ExactPoles& poles, const double* residues,
                             const Eigen::VectorXd& w);

/// sum_p a_p f(w - z_p) with given weights a_p (e.g. r_n1 r_n2)
Eigen::VectorXd ValuesDirectWeighted(const ExactPoles& poles,
                                     const double* weights,
                                     const Eigen::VectorXd& w);

/// Off-diagonal elements for all pairs of the Q levels whose residues are
/// the columns of residues_all (npoles x Q):
///   O(n1,n2) = sum_p r_n1(p) r_n2(p) [f(w_n1 - z_p) + f(w_n2 - z_p)],
/// f(t) = t / (t^2 + eta^2), diagonal zero, as blocked matrix products.
Eigen::MatrixXd OffDiagonal(const ExactPoles& poles,
                            const Eigen::MatrixXd& residues_all,
                            const Eigen::VectorXd& w);

}  // namespace exact_pole_sums
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EXACT_POLE_SUMS_H
