/*
 *            Copyright 2009-2023 The VOTCA Development Team
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
#ifndef VOTCA_XTP_QP_SOLVER_UTILS_H
#define VOTCA_XTP_QP_SOLVER_UTILS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/optional.hpp>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {
namespace qp_solver {

// Batched evaluation hook: a QP function may offer Prefetch(frequencies) to
// evaluate several points at once; others are evaluated point by point.
template <typename QPFunc>
auto PrefetchIfAvailable(const QPFunc& fqp, const std::vector<double>& nodes,
                         int) -> decltype(fqp.Prefetch(nodes), void()) {
  fqp.Prefetch(nodes);
}
template <typename QPFunc>
void PrefetchIfAvailable(const QPFunc&, const std::vector<double>&, long) {}

enum class EvalStage { Scan, Refine, Derivative, Other };

struct Stats {
  std::size_t sigma_scan_calls = 0;
  std::size_t sigma_refine_calls = 0;
  std::size_t sigma_derivative_calls = 0;
  std::size_t sigma_other_calls = 0;

  std::size_t sigma_repeat_calls = 0;
  std::size_t sigma_unique_frequencies = 0;

  std::size_t deriv_calls = 0;

  void Add(const Stats& other) {
    sigma_scan_calls += other.sigma_scan_calls;
    sigma_refine_calls += other.sigma_refine_calls;
    sigma_derivative_calls += other.sigma_derivative_calls;
    sigma_other_calls += other.sigma_other_calls;
    sigma_repeat_calls += other.sigma_repeat_calls;
    sigma_unique_frequencies += other.sigma_unique_frequencies;
    deriv_calls += other.deriv_calls;
  }

  std::size_t TotalSigmaCalls() const {
    return sigma_scan_calls + sigma_refine_calls + sigma_derivative_calls +
           sigma_other_calls;
  }
};

struct RootCandidate {
  double omega = 0.0;
  double residual = 0.0;
  double deriv = 0.0;
  double Z = 0.0;
  double distance_to_ref = 0.0;
  bool accepted = false;
};

struct WindowDiagnostics {
  Index shells_explored = 0;
  Index first_interval_shell = -1;
  Index first_accepted_shell = -1;
  Index chosen_shell = -1;
  Index intervals_found = 0;
};

struct SolverOptions {
  double g_sc_limit = 1e-5;
  Index qp_bisection_max_iter = 200;

  // The new grid-search design separates four previously entangled concerns:
  //   (1) full QP search extent
  //   (2) dense sign-change detection spacing
  //   (3) adaptive shell spacing / shell count
  //   (4) interval refinement method (configured elsewhere)
  // This makes the search numerics easier to reason about and allows robust
  // dense fallback without forcing the adaptive scan to inherit the same mesh.
  //
  // qp_full_window_half_width defines the full symmetric interval around the
  // current reference frequency that can be searched when no narrower window
  // is imposed by the caller.
  double qp_full_window_half_width = 0.75;

  // qp_dense_spacing is only for the robust dense scan that detects sign
  // changes before interval refinement. It is intentionally independent from
  // the adaptive search spacing.
  double qp_dense_spacing = 0.002;

  // The adaptive search starts at the reference energy and expands outward in
  // shells. The shell spacing can be set directly or derived from a requested
  // number of shells per half-window. This search is meant to find the most
  // relevant root early, while the dense scan remains the robust fallback.
  double qp_adaptive_shell_width = 0.025;
  Index qp_adaptive_shell_count = 0;  // 0 -> use width

  double min_accepted_Z = 0.05;
  double max_accepted_Z = 1.5;

  // Among accepted roots take the one closest to the reference frequency
  // (the previous iteration's solution) instead of the largest weight; set
  // for evGW iterations after the first with gw.qp_root_continuity.
  bool prefer_nearest_root = false;
};

// Legacy mapping helpers preserve the old qp_grid_steps / qp_grid_spacing
// behavior for existing XML files and tests. They are translation utilities
// only; the actual search code below works exclusively with the decoupled
// canonical options.
template <typename Opt>
inline double LegacyFullWindowHalfWidth(const Opt& opt) {
  if (opt.qp_grid_steps <= 1 || opt.qp_grid_spacing <= 0.0) {
    return -1.0;
  }
  return 0.5 * opt.qp_grid_spacing * double(opt.qp_grid_steps - 1);
}

template <typename Opt>
inline double LegacyAdaptiveShellWidth(const Opt& opt) {
  if (opt.qp_grid_steps <= 1 || opt.qp_grid_spacing <= 0.0) {
    return -1.0;
  }

  const double full_window_width =
      opt.qp_grid_spacing * double(opt.qp_grid_steps - 1);

  // Reconstruct the historical "coarse" spacing as closely as possible.
  // The hard-coded 21-point floor is legacy behavior from the original
  // implementation and intentionally remains confined to this compatibility
  // path instead of the new adaptive search itself.
  const Index base_coarse_steps = std::max<Index>(21, opt.qp_grid_steps / 4);

  if (base_coarse_steps <= 1) {
    return 4.0 * opt.qp_grid_spacing;
  }

  return full_window_width / double(base_coarse_steps - 1);
}

// NormalizeGridSearchOptions is the single entry point that converts a mix of
// new and legacy settings into one canonical set of search controls. New-style
// options always take precedence; legacy qp_grid_steps / qp_grid_spacing only
// fill values that were left unset by the caller.
template <typename Opt>
inline void NormalizeGridSearchOptions(Opt& opt) {
  const bool has_legacy = (opt.qp_grid_steps > 1 && opt.qp_grid_spacing > 0.0);

  if (opt.qp_full_window_half_width <= 0.0) {
    opt.qp_full_window_half_width =
        has_legacy ? LegacyFullWindowHalfWidth(opt) : 0.75;
  }

  if (opt.qp_dense_spacing <= 0.0) {
    opt.qp_dense_spacing = has_legacy ? opt.qp_grid_spacing : 0.002;
  }

  if (opt.qp_adaptive_shell_count <= 0 && opt.qp_adaptive_shell_width <= 0.0) {
    opt.qp_adaptive_shell_width =
        has_legacy ? LegacyAdaptiveShellWidth(opt) : 0.025;
  }

  if (opt.qp_full_window_half_width <= 0.0) {
    throw std::runtime_error(
        "Invalid QP search setup: qp_full_window_half_width must be > 0");
  }

  if (opt.qp_dense_spacing <= 0.0) {
    throw std::runtime_error(
        "Invalid QP search setup: qp_dense_spacing must be > 0");
  }

  if (opt.qp_adaptive_shell_count <= 0 && opt.qp_adaptive_shell_width <= 0.0) {
    throw std::runtime_error(
        "Invalid QP search setup: need qp_adaptive_shell_width > 0 or "
        "qp_adaptive_shell_count > 0");
  }
}

// Convert either an explicit shell width or an explicit shell count into the
// actual shell spacing used by the adaptive outward scan.
inline double EffectiveAdaptiveShellWidth(const SolverOptions& opt) {
  if (opt.qp_adaptive_shell_count > 0) {
    return opt.qp_full_window_half_width /
           static_cast<double>(opt.qp_adaptive_shell_count);
  }
  return opt.qp_adaptive_shell_width;
}

template <typename QPFunc>
double SolveQP_Bisection(double lowerbound, double f_lowerbound,
                         double upperbound, double f_upperbound,
                         const QPFunc& f, const SolverOptions& opt) {
  if (f_lowerbound * f_upperbound > 0.0) {
    throw std::runtime_error(
        "Bisection needs a positive and negative function value");
  }

  while (true) {
    const double c = 0.5 * (lowerbound + upperbound);
    if (std::abs(upperbound - lowerbound) < opt.g_sc_limit) {
      return c;
    }

    const double y_c = f.value(c, EvalStage::Refine);
    if (std::abs(y_c) < opt.g_sc_limit) {
      return c;
    }

    if (y_c * f_lowerbound > 0.0) {
      lowerbound = c;
      f_lowerbound = y_c;
    } else {
      upperbound = c;
      f_upperbound = y_c;
    }
  }
}

template <typename QPFunc>
double SolveQP_Brent(double lowerbound, double f_lowerbound, double upperbound,
                     double f_upperbound, const QPFunc& f,
                     const SolverOptions& opt) {
  if (f_lowerbound * f_upperbound > 0.0) {
    throw std::runtime_error(
        "Brent needs a positive and negative function value");
  }

  double a = lowerbound;
  double b = upperbound;
  double fa = f_lowerbound;
  double fb = f_upperbound;

  // c is the previous bracket endpoint
  double c = a;
  double fc = fa;

  // d and e track the last and second-last step sizes
  double d = b - a;
  double e = d;

  for (Index iter = 0; iter < opt.qp_bisection_max_iter; ++iter) {

    // Ensure that b is the best current estimate
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a;
      fc = fa;
      d = b - a;
      e = d;
    }

    if (std::abs(fc) < std::abs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    const double tol = opt.g_sc_limit;
    const double m = 0.5 * (c - b);

    if (std::abs(m) < tol || std::abs(fb) < opt.g_sc_limit) {
      return b;
    }

    if (std::abs(e) >= tol && std::abs(fa) > std::abs(fb)) {
      // Attempt inverse interpolation
      double s = fb / fa;
      double p = 0.0;
      double q = 0.0;

      if (a == c) {
        // Secant step
        p = 2.0 * m * s;
        q = 1.0 - s;
      } else {
        // Inverse quadratic interpolation
        double q1 = fa / fc;
        double r = fb / fc;
        p = s * (2.0 * m * q1 * (q1 - r) - (b - a) * (r - 1.0));
        q = (q1 - 1.0) * (r - 1.0) * (s - 1.0);
      }

      if (p > 0.0) {
        q = -q;
      }
      p = std::abs(p);

      // Accept interpolation only if it is safe
      if (q != 0.0 && 2.0 * p < std::min(3.0 * m * q - std::abs(tol * q),
                                         std::abs(e * q))) {
        e = d;
        d = p / q;
      } else {
        d = m;
        e = m;
      }
    } else {
      d = m;
      e = m;
    }

    a = b;
    fa = fb;

    if (std::abs(d) > tol) {
      b += d;
    } else {
      b += (m > 0.0 ? tol : -tol);
    }

    fb = f.value(b, EvalStage::Refine);
  }

  throw std::runtime_error(
      "Brent did not converge within qp_bisection_max_iter");
}

inline bool AcceptRoot(const RootCandidate& cand, const SolverOptions& opt) {
  if (!std::isfinite(cand.omega) || !std::isfinite(cand.Z)) {
    return false;
  }

  if (std::abs(cand.residual) > opt.g_sc_limit) {
    return false;
  }

  if (cand.Z <= 0.0) {
    return false;
  }

  if (cand.Z < opt.min_accepted_Z) {
    return false;
  }

  if (cand.Z > opt.max_accepted_Z) {
    return false;
  }

  return true;
}

inline double ScoreRoot(const RootCandidate& cand) {
  return cand.Z - 0.1 * cand.distance_to_ref;
}

/// The accepted root to use: the largest weight (ScoreRoot), or with
/// prefer_nearest_root the one closest to the reference. roots not empty.
inline const RootCandidate& SelectRoot(const std::vector<RootCandidate>& roots,
                                       const SolverOptions& opt) {
  if (opt.prefer_nearest_root) {
    return *std::min_element(
        roots.begin(), roots.end(),
        [](const RootCandidate& a, const RootCandidate& b) {
          return a.distance_to_ref < b.distance_to_ref;
        });
  }
  return *std::max_element(roots.begin(), roots.end(),
                           [](const RootCandidate& a, const RootCandidate& b) {
                             return ScoreRoot(a) < ScoreRoot(b);
                           });
}

/// Two accepted roots carry comparable weight when the second largest Z is
/// at least half the largest: the QP picture of the level is ambiguous and
/// the choice between them sensitive to small changes.
struct CompetingRoots {
  bool competing = false;
  double omega = 0.0;  // chosen root
  double Z = 0.0;
  double omega_alt = 0.0;  // the largest-weight other root
  double Z_alt = 0.0;
};

inline CompetingRoots CheckCompetingRoots(
    const std::vector<RootCandidate>& roots, double chosen) {
  CompetingRoots result;
  if (roots.size() < 2) {
    return result;
  }
  std::vector<RootCandidate> sorted = roots;
  std::sort(
      sorted.begin(), sorted.end(),
      [](const RootCandidate& a, const RootCandidate& b) { return a.Z > b.Z; });
  if (sorted[1].Z < 0.5 * sorted[0].Z) {
    return result;
  }
  result.competing = true;
  for (const RootCandidate& r : sorted) {
    if (r.omega == chosen) {
      result.omega = r.omega;
      result.Z = r.Z;
    } else if (result.Z_alt == 0.0) {
      result.omega_alt = r.omega;
      result.Z_alt = r.Z;
    }
  }
  return result;
}

template <typename QPFunc>
boost::optional<RootCandidate> RefineQPInterval(
    double lowerbound, double f_lowerbound, double upperbound,
    double f_upperbound, const QPFunc& f, double reference,
    const SolverOptions& opt, bool use_brent) {
  RootCandidate cand;

  // Near-exact-zero brackets: one endpoint is already essentially at the root.
  // Bisection and Brent both require opposite signs; handle this case directly
  // by returning whichever endpoint is closer to zero. This occurs when the
  // local mini-scan inside refine_and_store finds a near-zero node and inserts
  // a same-sign bracket -- those brackets are valid root indicators but cannot
  // be refined by interval methods.
  const bool left_near_zero = std::abs(f_lowerbound) <= opt.g_sc_limit;
  const bool right_near_zero = std::abs(f_upperbound) <= opt.g_sc_limit;
  const bool same_sign = (f_lowerbound * f_upperbound > 0.0);

  if (same_sign) {
    if (left_near_zero || right_near_zero) {
      // Pick the endpoint closest to zero as the root estimate
      cand.omega = (std::abs(f_lowerbound) <= std::abs(f_upperbound))
                       ? lowerbound
                       : upperbound;
    } else {
      // No bracket and no near-zero endpoint -- cannot refine, skip silently
      return boost::none;
    }
  } else {
    cand.omega = use_brent
                     ? SolveQP_Brent(lowerbound, f_lowerbound, upperbound,
                                     f_upperbound, f, opt)
                     : SolveQP_Bisection(lowerbound, f_lowerbound, upperbound,
                                         f_upperbound, f, opt);
  }

  cand.residual = f.value(cand.omega, EvalStage::Refine);
  cand.deriv = f.deriv(cand.omega);

  if (std::abs(cand.deriv) > 1e-14) {
    cand.Z = -1.0 / cand.deriv;
  } else {
    cand.Z = std::numeric_limits<double>::infinity();
  }

  cand.distance_to_ref = std::abs(cand.omega - reference);
  cand.accepted = AcceptRoot(cand, opt);

  cand.accepted = AcceptRoot(cand, opt);

  /*std::cout << "DEBUG ROOT "
            << "omega=" << cand.omega
            << " residual=" << cand.residual
            << " deriv=" << cand.deriv
            << " Z=" << cand.Z
            << " accepted=" << cand.accepted
            << std::endl;*/
  return cand;
}

/// Root tracking: the root of the branch found in the previous iteration,
/// searched only in [center - halfwidth, center + halfwidth] with the dense
/// spacing. Returns the accepted root closest to center whose weight is at
/// least min_Z_fraction of the previous one (the branch has not faded into a
/// satellite); none if there is no such root, so that the caller falls back
/// to the full search.
template <typename QPFunc>
boost::optional<RootCandidate> SolveQP_Tracked(
    QPFunc& fqp, double center, double halfwidth, double Z_previous,
    double min_Z_fraction, const SolverOptions& opt, bool use_brent = false) {
  const double spacing = opt.qp_dense_spacing;
  const Index n_steps = std::max<Index>(
      2, static_cast<Index>(std::ceil(2.0 * halfwidth / spacing)) + 1);
  const double left = center - halfwidth;
  const double dx = 2.0 * halfwidth / static_cast<double>(n_steps - 1);
  std::vector<double> nodes(static_cast<std::size_t>(n_steps));
  for (Index i = 0; i < n_steps; ++i) {
    nodes[std::size_t(i)] = left + static_cast<double>(i) * dx;
  }
  PrefetchIfAvailable(fqp, nodes, 0);

  boost::optional<RootCandidate> best;
  double f_prev = fqp.value(nodes[0], EvalStage::Scan);
  for (Index i = 1; i < n_steps; ++i) {
    const double x = nodes[std::size_t(i)];
    const double f = fqp.value(x, EvalStage::Scan);
    if (f_prev * f <= 0.0) {
      auto cand = RefineQPInterval(nodes[std::size_t(i - 1)], f_prev, x, f, fqp,
                                   center, opt, use_brent);
      if (cand && cand->accepted && cand->Z >= min_Z_fraction * Z_previous &&
          (!best || cand->distance_to_ref < best->distance_to_ref)) {
        best = cand;
      }
    }
    f_prev = f;
  }
  return best;
}

template <typename QPFunc>
boost::optional<double> SolveQP_Grid_Windowed(
    QPFunc& fqp, double frequency0, double left_limit, double right_limit,
    Index gw_sc_iteration, const SolverOptions& opt,
    WindowDiagnostics* wdiag = nullptr,
    std::vector<RootCandidate>* accepted_roots_out = nullptr,
    std::vector<RootCandidate>* rejected_roots_out = nullptr,
    bool use_brent = false) {
  struct SamplePoIndex {
    double omega = 0.0;
    double fval = 0.0;
  };

  WindowDiagnostics local_diag;
  std::vector<RootCandidate> accepted_roots;
  std::vector<RootCandidate> rejected_roots;

  if (left_limit >= right_limit) {
    if (wdiag != nullptr) {
      *wdiag = local_diag;
    }
    if (accepted_roots_out != nullptr) {
      *accepted_roots_out = accepted_roots;
    }
    if (rejected_roots_out != nullptr) {
      *rejected_roots_out = rejected_roots;
    }
    return boost::none;
  }

  const double shell_width = EffectiveAdaptiveShellWidth(opt);

  // The adaptive search is centered on the current reference frequency. On
  // the first GW iteration, a local linearized estimate is used when it stays
  // inside the requested window; later iterations keep the current frequency as
  // the center. From that center, the scan expands symmetrically in shells.
  double center = frequency0;

  if (gw_sc_iteration == 0) {
    const double f0 = fqp.value(frequency0, EvalStage::Other);
    const double df0 = fqp.deriv(frequency0);
    if (std::isfinite(f0) && std::isfinite(df0) && std::abs(df0) > 1e-6) {
      const double w_lin = frequency0 - f0 / df0;
      if (std::isfinite(w_lin) && w_lin >= left_limit && w_lin <= right_limit) {
        center = w_lin;
      }
    }
  }

  center = std::max(left_limit, std::min(right_limit, center));

  // Compute how many shells are needed to reach both limits from the chosen
  // center. Unlike the legacy implementation, this depends only on the
  // explicit adaptive shell control, not on the dense scan spacing.
  const double max_shell_reach =
      std::max(center - left_limit, right_limit - center);
  const Index n_shells =
      static_cast<Index>(std::ceil(max_shell_reach / shell_width));

  auto refine_and_store = [&](double a, double fa, double b, double fb,
                              Index shell_idx) {
    if (b < a) {
      std::swap(a, b);
      std::swap(fa, fb);
    }

    // Deterministic local mini-scan inside the coarse sign-change interval.
    // This restores much of the old dense-grid robustness while keeping the
    // fast coarse search globally.
    struct LocalBracket {
      double left = 0.0;
      double f_left = 0.0;
      double right = 0.0;
      double f_right = 0.0;
      double midpoint() const { return 0.5 * (left + right); }
    };

    const Index local_substeps = 12;
    std::vector<LocalBracket> local_brackets;
    local_brackets.reserve(static_cast<std::size_t>(local_substeps));

    if (b > a) {
      const double dx = (b - a) / static_cast<double>(local_substeps);
      {
        std::vector<double> nodes;
        for (Index i = 1; i < local_substeps; ++i) {
          nodes.push_back(a + static_cast<double>(i) * dx);
        }
        PrefetchIfAvailable(fqp, nodes, 0);
      }

      double x_prev = a;
      double f_prev = fa;

      for (Index i = 1; i <= local_substeps; ++i) {
        const double x_curr =
            (i == local_substeps) ? b : (a + static_cast<double>(i) * dx);
        const double f_curr =
            (i == local_substeps) ? fb : fqp.value(x_curr, EvalStage::Scan);

        // Standard sign change
        if ((f_prev < 0.0 && f_curr > 0.0) || (f_prev > 0.0 && f_curr < 0.0)) {
          local_brackets.push_back({x_prev, f_prev, x_curr, f_curr});
        }

        // Near-exact zero at the left node
        if (std::abs(f_prev) <= opt.g_sc_limit && x_prev < x_curr) {
          local_brackets.push_back({x_prev, f_prev, x_curr, f_curr});
        }

        // Near-exact zero at the right node
        if (std::abs(f_curr) <= opt.g_sc_limit && x_prev < x_curr) {
          local_brackets.push_back({x_prev, f_prev, x_curr, f_curr});
        }

        x_prev = x_curr;
        f_prev = f_curr;
      }
    }

    // Choose the local bracket deterministically:
    // midpoint closest to the search center, tie-broken by lower left edge.
    if (!local_brackets.empty()) {
      auto best_it = local_brackets.begin();
      double best_dist = std::abs(best_it->midpoint() - center);

      for (auto it = local_brackets.begin() + 1; it != local_brackets.end();
           ++it) {
        const double dist = std::abs(it->midpoint() - center);
        if (dist < best_dist - 1e-14 ||
            (std::abs(dist - best_dist) <= 1e-14 && it->left < best_it->left)) {
          best_it = it;
          best_dist = dist;
        }
      }

      a = best_it->left;
      fa = best_it->f_left;
      b = best_it->right;
      fb = best_it->f_right;
    }

    auto cand_opt =
        RefineQPInterval(a, fa, b, fb, fqp, frequency0, opt, use_brent);
    if (!cand_opt) {
      return;
    }

    if (local_diag.first_interval_shell < 0) {
      local_diag.first_interval_shell = shell_idx;
    }
    ++local_diag.intervals_found;

    const RootCandidate& cand = cand_opt.value();
    if (cand.accepted) {
      if (local_diag.first_accepted_shell < 0) {
        local_diag.first_accepted_shell = shell_idx;
      }
      accepted_roots.push_back(cand);
    } else {
      rejected_roots.push_back(cand);
    }
  };

  // All shell points are known in advance: evaluate them in one batch.
  {
    std::vector<double> nodes{center};
    for (Index shell = 1; shell <= n_shells; ++shell) {
      const double delta = double(shell) * shell_width;
      if (center - delta >= left_limit) {
        nodes.push_back(center - delta);
      }
      if (center + delta <= right_limit) {
        nodes.push_back(center + delta);
      }
    }
    nodes.push_back(left_limit);
    nodes.push_back(right_limit);
    PrefetchIfAvailable(fqp, nodes, 0);
  }

  SamplePoIndex center_pt{center, fqp.value(center, EvalStage::Scan)};

  bool left_active = true;
  bool right_active = true;
  SamplePoIndex left_prev = center_pt;
  SamplePoIndex right_prev = center_pt;

  for (Index shell = 1; shell <= n_shells; ++shell) {
    local_diag.shells_explored = shell;
    bool added_this_shell = false;
    const double delta = double(shell) * shell_width;

    // Each shell adds one point to the left and one to the right of the
    // center. Whenever a sign change is detected between neighboring samples,
    // that interval is refined by bisection or Brent and the resulting root is
    // classified as accepted or rejected using the Z-based acceptance test.

    if (left_active) {
      const double omega_left = center - delta;
      if (omega_left >= left_limit) {
        SamplePoIndex left_curr{omega_left,
                                fqp.value(omega_left, EvalStage::Scan)};
        added_this_shell = true;

        if (left_prev.fval * left_curr.fval < 0.0) {
          refine_and_store(left_curr.omega, left_curr.fval, left_prev.omega,
                           left_prev.fval, shell);
        }
        left_prev = left_curr;
      } else {
        left_active = false;
      }
    }

    if (right_active) {
      const double omega_right = center + delta;
      if (omega_right <= right_limit) {
        SamplePoIndex right_curr{omega_right,
                                 fqp.value(omega_right, EvalStage::Scan)};
        added_this_shell = true;

        if (right_prev.fval * right_curr.fval < 0.0) {
          refine_and_store(right_prev.omega, right_prev.fval, right_curr.omega,
                           right_curr.fval, shell);
        }
        right_prev = right_curr;
      } else {
        right_active = false;
      }
    }

    if (!added_this_shell && !left_active && !right_active) {
      break;
    }
  }

  if (left_prev.omega > left_limit + 1e-12) {
    SamplePoIndex left_end{left_limit, fqp.value(left_limit, EvalStage::Scan)};
    if (left_end.fval * left_prev.fval < 0.0) {
      refine_and_store(left_end.omega, left_end.fval, left_prev.omega,
                       left_prev.fval, local_diag.shells_explored + 1);
    }
  }

  if (right_prev.omega < right_limit - 1e-12) {
    SamplePoIndex right_end{right_limit,
                            fqp.value(right_limit, EvalStage::Scan)};
    if (right_prev.fval * right_end.fval < 0.0) {
      refine_and_store(right_prev.omega, right_prev.fval, right_end.omega,
                       right_end.fval, local_diag.shells_explored + 1);
    }
  }

  // Preference order for the adaptive search:
  //   1. return the best accepted root found in the scanned shells
  //   2. if the caller allows it, return the best rejected root
  // The outer GW driver uses this to enforce that restricted-window searches
  // only succeed on accepted roots; otherwise it escalates to a full dense
  // scan.
  if (!accepted_roots.empty()) {
    const RootCandidate* best = &SelectRoot(accepted_roots, opt);

    local_diag.chosen_shell = static_cast<int>(
        std::llround(std::abs(best->omega - center) / shell_width));

    if (wdiag != nullptr) {
      *wdiag = local_diag;
    }
    if (accepted_roots_out != nullptr) {
      *accepted_roots_out = accepted_roots;
    }
    if (rejected_roots_out != nullptr) {
      *rejected_roots_out = rejected_roots;
    }
    return best->omega;
  }

  if (!rejected_roots.empty()) {
    auto least_bad =
        std::max_element(rejected_roots.begin(), rejected_roots.end(),
                         [](const RootCandidate& a, const RootCandidate& b) {
                           return ScoreRoot(a) < ScoreRoot(b);
                         });

    local_diag.chosen_shell = static_cast<int>(
        std::llround(std::abs(least_bad->omega - center) / shell_width));

    if (wdiag != nullptr) {
      *wdiag = local_diag;
    }
    if (accepted_roots_out != nullptr) {
      *accepted_roots_out = accepted_roots;
    }
    if (rejected_roots_out != nullptr) {
      *rejected_roots_out = rejected_roots;
    }
    return least_bad->omega;
  }

  if (wdiag != nullptr) {
    *wdiag = local_diag;
  }
  if (accepted_roots_out != nullptr) {
    *accepted_roots_out = accepted_roots;
  }
  if (rejected_roots_out != nullptr) {
    *rejected_roots_out = rejected_roots;
  }
  return boost::none;
}

/// The n levels with the largest QP residual |solution - input| (largest
/// first) as "level: input -> solution (residual)", for the evGW log.
/// Levels are numbered from first_level.
inline std::string LargestResiduals(const Eigen::VectorXd& solution,
                                    const Eigen::VectorXd& input,
                                    Index first_level, Index n = 5) {
  const Index size = std::min(solution.size(), input.size());
  std::vector<Index> order(static_cast<std::size_t>(size));
  for (Index i = 0; i < size; ++i) {
    order[std::size_t(i)] = i;
  }
  const Index shown = std::min(n, size);
  std::partial_sort(order.begin(), order.begin() + shown, order.end(),
                    [&](Index a, Index b) {
                      return std::abs(solution(a) - input(a)) >
                             std::abs(solution(b) - input(b));
                    });
  std::string out;
  char buffer[96];
  for (Index k = 0; k < shown; ++k) {
    const Index i = order[std::size_t(k)];
    std::snprintf(buffer, sizeof(buffer), "%s%ld: %+.6f -> %+.6f (%+.1e)",
                  (k == 0) ? "" : ", ", static_cast<long>(i + first_level),
                  input(i), solution(i), solution(i) - input(i));
    out += buffer;
  }
  return out;
}

}  // namespace qp_solver
}  // namespace xtp
}  // namespace votca

#endif