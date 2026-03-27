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
#include <limits>
#include <stdexcept>
#include <vector>

#include <boost/optional.hpp>

namespace votca {
namespace xtp {
namespace qp_solver {

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
  double qp_grid_spacing = 0.01;
  Index qp_grid_steps = 201;

  double min_accepted_Z = 0.05;
  double max_accepted_Z = 1.5;
};

template <typename QPFunc>
double SolveQP_Bisection(double lowerbound, double f_lowerbound,
                         double upperbound, double f_upperbound,
                         const QPFunc& f, const SolverOptions& opt) {
  if (f_lowerbound * f_upperbound > 0.0) {
    throw std::runtime_error(
        "Bisection needs a positive and negative function value");
  }

  std::cout << "Bisection" << std::endl;

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

  std::cout << "Brent" << std::endl;

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

template <typename QPFunc>
boost::optional<RootCandidate> RefineQPInterval(
    double lowerbound, double f_lowerbound, double upperbound,
    double f_upperbound, const QPFunc& f, double reference,
    const SolverOptions& opt, bool use_brent) {
  RootCandidate cand;

  cand.omega = use_brent ? SolveQP_Brent(lowerbound, f_lowerbound, upperbound,
                                         f_upperbound, f, opt)
                         : SolveQP_Bisection(lowerbound, f_lowerbound,
                                             upperbound, f_upperbound, f, opt);

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

  const double range = right_limit - left_limit;
  const double full_window_width =
      opt.qp_grid_spacing * double(opt.qp_grid_steps - 1);
  const Index base_coarse_steps = std::max<Index>(21, opt.qp_grid_steps / 4);
  const double base_coarse_spacing =
      full_window_width / double(base_coarse_steps - 1);

  const Index coarse_steps =
      std::max(3, static_cast<int>(std::ceil(range / base_coarse_spacing)) + 1);
  const double coarse_spacing = range / double(coarse_steps - 1);

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

  auto refine_and_store = [&](double a, double fa, double b, double fb,
                              Index shell_idx) {
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

  SamplePoIndex center_pt{center, fqp.value(center, EvalStage::Scan)};

  bool left_active = true;
  bool right_active = true;
  SamplePoIndex left_prev = center_pt;
  SamplePoIndex right_prev = center_pt;

  for (Index shell = 1; shell < coarse_steps; ++shell) {
    local_diag.shells_explored = shell;
    bool added_this_shell = false;
    const double delta = double(shell) * coarse_spacing;

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

  if (!accepted_roots.empty()) {
    auto best =
        std::max_element(accepted_roots.begin(), accepted_roots.end(),
                         [](const RootCandidate& a, const RootCandidate& b) {
                           return ScoreRoot(a) < ScoreRoot(b);
                         });

    local_diag.chosen_shell = static_cast<int>(
        std::llround(std::abs(best->omega - center) / coarse_spacing));

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
        std::llround(std::abs(least_bad->omega - center) / coarse_spacing));

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

}  // namespace qp_solver
}  // namespace xtp
}  // namespace votca

#endif