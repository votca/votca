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
#ifndef VOTCA_XTP_OUT_OF_WINDOW_SHIFT_H
#define VOTCA_XTP_OUT_OF_WINDOW_SHIFT_H

// Standard includes
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief How the RPA energies of levels outside the GW window are shifted.
 *
 * Levels in the RPA but not in the GW window [qpmin, qpmax] get no QP
 * correction of their own. With D_i = E_i - e_i the corrections of the
 * window levels, occupied (below qpmin) and virtual (above qpmax) levels are
 * shifted by
 *
 *  - max:     the largest |D_i| of the occupied (virtual) window levels,
 *             downwards (upwards), one rigid shift per side;
 *  - nearest: D at the window boundary (qpmin, qpmax), with its sign;
 *  - linear:  D(e) = a + b (e - e_b) from a least-squares fit of D_i
 *             against e_i over the fraction fit_fraction of the occupied
 *             (virtual) window levels next to the boundary e_b; it is
 *             extrapolated at most as far beyond the boundary as the fitted
 *             levels span, and is constant beyond that.
 *
 * nearest and linear are smooth functions of the QP energies; max is not
 * (which level holds the maximum changes), which matters for the fixed-point
 * iterations of evGW and QSGW.
 */
struct OutOfWindowShift {
  enum class Mode { Max, Nearest, Linear };
  Mode mode = Mode::Max;
  double fit_fraction = 0.3;

  static Mode Parse(const std::string& name) {
    if (name == "max") {
      return Mode::Max;
    } else if (name == "nearest") {
      return Mode::Nearest;
    } else if (name == "linear") {
      return Mode::Linear;
    }
    throw std::runtime_error("Unknown out_of_window_shift '" + name +
                             "'; choices are max, nearest, linear");
  }
  static std::string Name(Mode mode) {
    switch (mode) {
      case Mode::Nearest:
        return "nearest";
      case Mode::Linear:
        return "linear";
      default:
        return "max";
    }
  }

  /// Shifts the levels outside [qpmin, qpmax] of energies (indexed from
  /// rpamin, window levels already QP energies); dft and the level numbers
  /// homo, qpmin, qpmax are absolute.
  void Apply(Eigen::VectorXd& energies, const Eigen::VectorXd& dft,
             Index rpamin, Index rpamax, Index homo, Index qpmin,
             Index qpmax) const {
    const Index lumo = homo + 1;
    if (qpmin > rpamin) {
      ShiftSide(energies, dft, rpamin, rpamin, qpmin - 1, qpmin,
                std::min(homo, qpmax), true);
    }
    if (qpmax < rpamax) {
      ShiftSide(energies, dft, rpamin, qpmax + 1, rpamax, std::max(lumo, qpmin),
                qpmax, false);
    }
  }

 private:
  // Shifts levels [out_first, out_last] using the corrections of the window
  // levels [win_first, win_last] of the same kind; below is true for the
  // occupied side (boundary win_first), false for the virtual side (boundary
  // win_last).
  void ShiftSide(Eigen::VectorXd& energies, const Eigen::VectorXd& dft,
                 Index rpamin, Index out_first, Index out_last, Index win_first,
                 Index win_last, bool below) const {
    const Index nwin = win_last - win_first + 1;
    auto correction = [&](Index level) {
      return energies(level - rpamin) - dft(level);
    };
    auto out = energies.segment(out_first - rpamin, out_last - out_first + 1);
    if (nwin <= 0) {
      return;
    }
    if (mode == Mode::Max) {
      double max_correction = 0.0;
      for (Index l = win_first; l <= win_last; ++l) {
        max_correction = std::max(max_correction, std::abs(correction(l)));
      }
      out.array() += below ? -max_correction : max_correction;
      return;
    }
    const Index boundary = below ? win_first : win_last;
    const Index nfit =
        (mode == Mode::Nearest)
            ? 1
            : std::min(nwin,
                       std::max<Index>(
                           2, Index(std::ceil(fit_fraction * double(nwin)))));
    if (nfit < 2) {
      out.array() += correction(boundary);
      return;
    }
    // least squares D = a + b x over the nfit levels next to the boundary,
    // x = e - e_boundary
    const Index fit_first = below ? win_first : win_last - nfit + 1;
    const double e_b = dft(boundary);
    double sx = 0, sy = 0, sxx = 0, sxy = 0, span = 0;
    for (Index l = fit_first; l < fit_first + nfit; ++l) {
      const double x = dft(l) - e_b;
      const double y = correction(l);
      sx += x;
      sy += y;
      sxx += x * x;
      sxy += x * y;
      span = std::max(span, std::abs(x));
    }
    const double n = double(nfit);
    const double var = sxx - sx * sx / n;
    const double b = (var > 1e-12 * n) ? (sxy - sx * sy / n) / var : 0.0;
    const double a = (sy - b * sx) / n;
    for (Index l = out_first; l <= out_last; ++l) {
      const double x = std::clamp(dft(l) - e_b, -span, span);
      energies(l - rpamin) += a + b * x;
    }
  }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_OUT_OF_WINDOW_SHIFT_H
