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
#ifndef VOTCA_XTP_SCALED_ERFC_H
#define VOTCA_XTP_SCALED_ERFC_H

// Standard includes
#include <cmath>

namespace votca {
namespace xtp {

/// exp(x^2) erfc(x) for x >= 0. The direct product overflows to inf * 0
/// (NaN) beyond x ~ 26.6; there the asymptotic series is used (relative
/// error below 1e-11 from x = 26 on).
inline double ScaledErfc(double x) {
  if (x < 26.0) {
    return std::exp(x * x) * std::erfc(x);
  }
  const double inv2 = 1.0 / (x * x);
  return (1.0 - 0.5 * inv2 + 0.75 * inv2 * inv2 - 1.875 * inv2 * inv2 * inv2) /
         (x * 1.7724538509055160273);  // x sqrt(pi)
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SCALED_ERFC_H
