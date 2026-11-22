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

// Internal header of libxtp (libint2_calls.cc, libint2_derivative_calls.cc):
// include it after the libint2 headers.

#pragma once
#ifndef VOTCA_XTP_GAUSSIANSITESHELLS_H
#define VOTCA_XTP_GAUSSIANSITESHELLS_H

#include <array>

#include "votca/xtp/eigen.h"

namespace votca {
namespace xtp {

// libint2 shells of a Gaussian charge (s), dipole (Cartesian p) and
// quadrupole (Cartesian d) of exponent beta = 1/width^2 at pos, and the
// factors that turn the L2-normalized primitives libint2 makes of them into
// the unit charge g = (beta/pi)^3/2 exp(-beta r^2), the unit dipole density
// 2 beta (r - R) g = -grad g and the quadrupole density
// c sum_ij Theta_ij d_i d_j g (c = kQuadrupoleScale) of AOGaussianSites:
//   s: N_s = (2 beta/pi)^3/4, charge N_s (pi/beta)^3/2;
//   p: N_p = N_s 2 sqrt(beta);  d: N_d = N_s 4 beta/sqrt(3) for all six
//   components (libint2 normalizes them like x^2).
struct GaussianSiteShells {
  // potential of the point quadrupole: c sum_ij Theta_ij d_i d_j (1/r)
  // = 3c sum_ij Theta_ij x_i x_j / r^5, which with c = 1/4 is the potential
  // of AOMultipole for Theta = StaticSite::CalculateCartesianMultipole()
  static constexpr double kQuadrupoleScale = 0.25;

  libint2::Shell s;
  libint2::Shell p;
  libint2::Shell d;
  double s_factor;
  double p_factor;
  double d_factor;
  GaussianSiteShells(const Eigen::Vector3d& pos, double width);
  // weights of the six Cartesian d components (xx, xy, xz, yy, yz, zz) for
  // sum_ij Theta_ij x_i x_j
  static std::array<double, 6> QuadrupoleWeights(const Eigen::Matrix3d& Q);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GAUSSIANSITESHELLS_H
