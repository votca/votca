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

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace votca {
namespace xtp {

/**
 * @brief Build an initial guess for full-BSE Davidson diagonalization
 *
 * Strategy:
 *  - rank excitation channels by local diagonal full-BSE estimate
 *      omega_i = sqrt(max(0, A_ii^2 - B_ii^2))
 *  - construct X-only basis vectors for robustness with HAM Davidson
 *
 * This avoids:
 *  - expensive TDA pre-diagonalization
 *  - unstable X/Y initial guesses
 *
 * Works for:
 *  - UKS full-BSE
 *  - RKS full-BSE (future reuse)
 */
inline Eigen::MatrixXd BuildFullBSEXRankedInitialGuess(
    const Eigen::VectorXd& adiag, const Eigen::VectorXd& bdiag, Index nroots) {

  if (adiag.size() != bdiag.size()) {
    throw std::runtime_error("BuildFullBSEXRankedInitialGuess: size mismatch.");
  }

  const Index n = adiag.size();
  const Index nguess = std::min<Index>(n, std::max<Index>(4 * nroots, 8));

  struct RankedMode {
    double omega;
    double a;
    Index idx;
  };

  std::vector<RankedMode> ranked;
  ranked.reserve(n);

  for (Index i = 0; i < n; ++i) {
    const double a = adiag(i);
    const double b = bdiag(i);

    const double disc = std::max(0.0, (a - b) * (a + b));
    const double omega = std::sqrt(disc);

    ranked.push_back({omega, a, i});
  }

  std::sort(ranked.begin(), ranked.end(),
            [](const RankedMode& lhs, const RankedMode& rhs) {
              if (lhs.omega != rhs.omega) {
                return lhs.omega < rhs.omega;
              }
              return lhs.a < rhs.a;
            });

  Eigen::MatrixXd guess = Eigen::MatrixXd::Zero(2 * n, nguess);

  // X-only basis vectors
  for (Index col = 0; col < nguess; ++col) {
    const Index i = ranked[col].idx;
    guess(i, col) = 1.0;
  }

  return guess;
}

}  // namespace xtp
}  // namespace votca