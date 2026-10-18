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
#ifndef VOTCA_XTP_GAUSSIANKERNEL_H
#define VOTCA_XTP_GAUSSIANKERNEL_H

// Standard includes
#include <cmath>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief Interaction of two normalized Gaussian charge distributions
 * exp(-r^2/s_i^2) and exp(-r^2/s_j^2): erf(r/S)/r with S = sqrt(s_i^2 +
 * s_j^2), and its radial derivative functions
 *
 *   B0 = erf(r/S)/r
 *   B1 = -(1/r) dB0/dr = (B0 - 2/(sqrt(pi) S) exp(-r^2/S^2)) / r^2
 *   B2 = -(1/r) dB1/dr = (3 B1 - 4/(sqrt(pi) S^3) exp(-r^2/S^2)) / r^2
 *
 * S = 0 gives the bare Coulomb functions 1/r, 1/r^3, 3/r^5. With r the
 * vector from the source to the receiving site:
 *   potential of a charge q:   q B0      field of a charge q:   q B1 r
 *   potential of a dipole mu:  B1 mu.r   field of a dipole mu:  (B2 r r^T -
 *   B1 1) mu
 * (the Gaussian dipole is the derivative of a Gaussian charge). Used by the
 * charge-response model (charge-response sites of width s, Thole sites as
 * points).
 */
struct GaussianKernel {
  double B0 = 0.0;
  double B1 = 0.0;
  double B2 = 0.0;

  GaussianKernel(double r, double S) {
    if (S <= 0.0) {
      B0 = 1.0 / r;
      B1 = B0 * B0 * B0;
      B2 = 3.0 * B1 * B0 * B0;
      return;
    }
    const double x = r / S;
    constexpr double two_over_sqrtpi = 1.1283791670955126;
    if (x < 0.05) {
      // series in x (to x^8): avoids the cancellation of the closed form at
      // short range
      const double x2 = x * x;
      const double s1 = two_over_sqrtpi / S;
      B0 = s1 *
           (1.0 + x2 * (-1.0 / 3.0 +
                        x2 * (1.0 / 10.0 + x2 * (-1.0 / 42.0 + x2 / 216.0))));
      const double s3 = s1 / (S * S);
      B1 = s3 * (2.0 / 3.0 +
                 x2 * (-2.0 / 5.0 +
                       x2 * (1.0 / 7.0 + x2 * (-1.0 / 27.0 + x2 / 132.0))));
      const double s5 = s3 / (S * S);
      B2 = s5 * (4.0 / 5.0 +
                 x2 * (-4.0 / 7.0 +
                       x2 * (2.0 / 9.0 + x2 * (-2.0 / 33.0 + x2 / 78.0))));
      return;
    }
    const double g = std::exp(-x * x);
    B0 = std::erf(x) / r;
    const double r2 = r * r;
    B1 = (B0 - two_over_sqrtpi * g / S) / r2;
    B2 = (3.0 * B1 - 2.0 * two_over_sqrtpi * g / (S * S * S)) / r2;
  }

  /// gradient of the potential (= -field) tensor of a dipole at the source:
  /// grad phi = T mu with T = B1 1 - B2 r r^T
  Eigen::Matrix3d DipoleTensor(const Eigen::Vector3d& r) const {
    Eigen::Matrix3d T = -B2 * r * r.transpose();
    T.diagonal().array() += B1;
    return T;
  }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GAUSSIANKERNEL_H
