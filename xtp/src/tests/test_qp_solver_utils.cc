/*
 * Copyright 2009-2023 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE qp_solver_utils_test

#include <cmath>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "votca/tools/types.h"
#include "votca/xtp/qp_solver_utils.h"

using votca::Index;
using namespace votca::xtp::qp_solver;

namespace {

constexpr double kTol = 1e-10;

struct DummyOpt {
  Index qp_grid_steps = 0;
  double qp_grid_spacing = 0.0;

  double qp_full_window_half_width = -1.0;
  double qp_dense_spacing = -1.0;
  double qp_adaptive_shell_width = -1.0;
  Index qp_adaptive_shell_count = 0;
};

struct LinearAcceptedFunc {
  explicit LinearAcceptedFunc(double root_) : root(root_) {}

  double value(double omega, EvalStage /*stage*/) const { return root - omega; }

  double deriv(double /*omega*/) const { return -1.0; }

  double root;
};

struct LinearRejectedFunc {
  explicit LinearRejectedFunc(double root_) : root(root_) {}

  double value(double omega, EvalStage /*stage*/) const { return omega - root; }

  double deriv(double /*omega*/) const { return 1.0; }

  double root;
};

struct TwoRootAcceptedFunc {
  TwoRootAcceptedFunc(double root1_, double root2_)
      : root1(root1_), root2(root2_) {}

  double value(double omega, EvalStage /*stage*/) const {
    return (omega - root1) * (omega - root2);
  }

  // Deliberately constant and negative so both refined roots are "accepted"
  // by the Z criterion. This is a mock unit-test function meant to isolate
  // the root-selection logic, not a physical self-energy model.
  double deriv(double /*omega*/) const { return -1.0; }

  double root1;
  double root2;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(qp_solver_utils_test)

BOOST_AUTO_TEST_CASE(legacy_full_window_half_width_mapping) {
  DummyOpt opt;
  opt.qp_grid_steps = 201;
  opt.qp_grid_spacing = 0.01;

  const double half_width = LegacyFullWindowHalfWidth(opt);
  BOOST_CHECK_CLOSE_FRACTION(half_width, 1.0, kTol);
}

BOOST_AUTO_TEST_CASE(legacy_adaptive_shell_width_mapping) {
  DummyOpt opt;
  opt.qp_grid_steps = 201;
  opt.qp_grid_spacing = 0.01;

  // full width = (201 - 1) * 0.01 = 2.0
  // base_coarse_steps = max(21, 201 / 4) = max(21, 50) = 50
  // legacy shell width = 2.0 / (50 - 1)
  const double shell_width = LegacyAdaptiveShellWidth(opt);
  BOOST_CHECK_CLOSE_FRACTION(shell_width, 2.0 / 49.0, kTol);
}

BOOST_AUTO_TEST_CASE(normalize_uses_new_defaults_when_nothing_is_set) {
  DummyOpt opt;

  NormalizeGridSearchOptions(opt);

  BOOST_CHECK_CLOSE_FRACTION(opt.qp_full_window_half_width, 0.75, kTol);
  BOOST_CHECK_CLOSE_FRACTION(opt.qp_dense_spacing, 0.002, kTol);
  BOOST_CHECK_CLOSE_FRACTION(opt.qp_adaptive_shell_width, 0.025, kTol);
  BOOST_CHECK_EQUAL(opt.qp_adaptive_shell_count, 0);
}

BOOST_AUTO_TEST_CASE(normalize_legacy_values_fill_unset_new_controls) {
  DummyOpt opt;
  opt.qp_grid_steps = 201;
  opt.qp_grid_spacing = 0.01;

  NormalizeGridSearchOptions(opt);

  BOOST_CHECK_CLOSE_FRACTION(opt.qp_full_window_half_width, 1.0, kTol);
  BOOST_CHECK_CLOSE_FRACTION(opt.qp_dense_spacing, 0.01, kTol);
  BOOST_CHECK_CLOSE_FRACTION(opt.qp_adaptive_shell_width, 2.0 / 49.0, kTol);
  BOOST_CHECK_EQUAL(opt.qp_adaptive_shell_count, 0);
}

BOOST_AUTO_TEST_CASE(normalize_preserves_explicit_new_values_over_legacy) {
  DummyOpt opt;
  opt.qp_grid_steps = 201;
  opt.qp_grid_spacing = 0.01;

  opt.qp_full_window_half_width = 0.75;
  opt.qp_dense_spacing = 0.002;
  opt.qp_adaptive_shell_width = 0.03;
  opt.qp_adaptive_shell_count = 0;

  NormalizeGridSearchOptions(opt);

  BOOST_CHECK_CLOSE_FRACTION(opt.qp_full_window_half_width, 0.75, kTol);
  BOOST_CHECK_CLOSE_FRACTION(opt.qp_dense_spacing, 0.002, kTol);
  BOOST_CHECK_CLOSE_FRACTION(opt.qp_adaptive_shell_width, 0.03, kTol);
  BOOST_CHECK_EQUAL(opt.qp_adaptive_shell_count, 0);
}

BOOST_AUTO_TEST_CASE(shell_count_overrides_shell_width) {
  SolverOptions opt;
  opt.qp_full_window_half_width = 0.75;
  opt.qp_adaptive_shell_width = 0.03;
  opt.qp_adaptive_shell_count = 30;

  const double shell_width = EffectiveAdaptiveShellWidth(opt);
  BOOST_CHECK_CLOSE_FRACTION(shell_width, 0.75 / 30.0, kTol);
}

BOOST_AUTO_TEST_CASE(accept_root_applies_residual_and_Z_criteria) {
  SolverOptions opt;
  opt.g_sc_limit = 1e-5;
  opt.min_accepted_Z = 0.05;
  opt.max_accepted_Z = 1.5;

  RootCandidate good;
  good.omega = 0.12;
  good.residual = 1e-7;
  good.deriv = -1.25;
  good.Z = 0.8;
  good.distance_to_ref = 0.12;
  BOOST_CHECK(AcceptRoot(good, opt));

  RootCandidate bad_residual = good;
  bad_residual.residual = 1e-3;
  BOOST_CHECK(!AcceptRoot(bad_residual, opt));

  RootCandidate bad_small_Z = good;
  bad_small_Z.Z = 0.01;
  BOOST_CHECK(!AcceptRoot(bad_small_Z, opt));

  RootCandidate bad_large_Z = good;
  bad_large_Z.Z = 2.0;
  BOOST_CHECK(!AcceptRoot(bad_large_Z, opt));

  RootCandidate bad_negative_Z = good;
  bad_negative_Z.Z = -0.5;
  BOOST_CHECK(!AcceptRoot(bad_negative_Z, opt));
}

BOOST_AUTO_TEST_CASE(windowed_solver_finds_simple_root_with_bisection) {
  LinearAcceptedFunc fqp(0.23);

  SolverOptions opt;
  opt.g_sc_limit = 1e-8;
  opt.qp_bisection_max_iter = 200;
  opt.qp_full_window_half_width = 0.5;
  opt.qp_dense_spacing = 0.002;
  opt.qp_adaptive_shell_width = 0.05;
  opt.qp_adaptive_shell_count = 0;

  WindowDiagnostics diag;
  std::vector<RootCandidate> accepted_roots;
  std::vector<RootCandidate> rejected_roots;

  boost::optional<double> root =
      SolveQP_Grid_Windowed(fqp, 0.0, -0.5, 0.5, 1, opt, &diag, &accepted_roots,
                            &rejected_roots, false);

  BOOST_REQUIRE(root);
  BOOST_CHECK_CLOSE_FRACTION(*root, 0.23, 1e-6);

  BOOST_CHECK_EQUAL(accepted_roots.size(), 1);
  BOOST_CHECK_EQUAL(rejected_roots.size(), 0);

  BOOST_CHECK(accepted_roots.front().accepted);
  BOOST_CHECK_CLOSE_FRACTION(accepted_roots.front().Z, 1.0, kTol);

  BOOST_CHECK_EQUAL(diag.first_interval_shell, 5);
  BOOST_CHECK_EQUAL(diag.first_accepted_shell, 5);
  BOOST_CHECK_EQUAL(diag.chosen_shell, 5);
  BOOST_CHECK_EQUAL(diag.intervals_found, 1);
  BOOST_CHECK_GE(diag.shells_explored, 5);
}

BOOST_AUTO_TEST_CASE(windowed_solver_finds_simple_root_with_brent) {
  LinearAcceptedFunc fqp(0.23);

  SolverOptions opt;
  opt.g_sc_limit = 1e-8;
  opt.qp_bisection_max_iter = 200;
  opt.qp_full_window_half_width = 0.5;
  opt.qp_dense_spacing = 0.002;
  opt.qp_adaptive_shell_width = 0.05;
  opt.qp_adaptive_shell_count = 0;

  WindowDiagnostics diag;
  std::vector<RootCandidate> accepted_roots;
  std::vector<RootCandidate> rejected_roots;

  boost::optional<double> root =
      SolveQP_Grid_Windowed(fqp, 0.0, -0.5, 0.5, 1, opt, &diag, &accepted_roots,
                            &rejected_roots, true);

  BOOST_REQUIRE(root);
  BOOST_CHECK_CLOSE_FRACTION(*root, 0.23, 1e-10);

  BOOST_CHECK_EQUAL(accepted_roots.size(), 1);
  BOOST_CHECK_EQUAL(rejected_roots.size(), 0);
  BOOST_CHECK_EQUAL(diag.intervals_found, 1);
  BOOST_CHECK_GE(diag.shells_explored, 5);
}

BOOST_AUTO_TEST_CASE(windowed_solver_returns_nearest_accepted_root) {
  TwoRootAcceptedFunc fqp(0.12, 0.62);

  SolverOptions opt;
  opt.g_sc_limit = 1e-8;
  opt.qp_bisection_max_iter = 200;
  opt.qp_full_window_half_width = 0.8;
  opt.qp_dense_spacing = 0.002;
  opt.qp_adaptive_shell_width = 0.05;
  opt.qp_adaptive_shell_count = 0;

  WindowDiagnostics diag;
  std::vector<RootCandidate> accepted_roots;
  std::vector<RootCandidate> rejected_roots;

  boost::optional<double> root =
      SolveQP_Grid_Windowed(fqp, 0.0, -0.2, 0.8, 1, opt, &diag, &accepted_roots,
                            &rejected_roots, false);

  BOOST_REQUIRE(root);
  BOOST_CHECK_CLOSE_FRACTION(*root, 0.12, 1e-6);

  BOOST_CHECK_EQUAL(accepted_roots.size(), 2);
  BOOST_CHECK_EQUAL(rejected_roots.size(), 0);
  BOOST_CHECK_GE(diag.intervals_found, 2);
  BOOST_CHECK(diag.first_accepted_shell >= 0);
  BOOST_CHECK(diag.chosen_shell >= 0);
}

BOOST_AUTO_TEST_CASE(
    windowed_solver_returns_rejected_root_if_no_accepted_root_exists) {
  LinearRejectedFunc fqp(0.23);

  SolverOptions opt;
  opt.g_sc_limit = 1e-8;
  opt.qp_bisection_max_iter = 200;
  opt.qp_full_window_half_width = 0.5;
  opt.qp_dense_spacing = 0.002;
  opt.qp_adaptive_shell_width = 0.05;
  opt.qp_adaptive_shell_count = 0;

  WindowDiagnostics diag;
  std::vector<RootCandidate> accepted_roots;
  std::vector<RootCandidate> rejected_roots;

  boost::optional<double> root =
      SolveQP_Grid_Windowed(fqp, 0.0, -0.5, 0.5, 1, opt, &diag, &accepted_roots,
                            &rejected_roots, false);

  BOOST_REQUIRE(root);
  BOOST_CHECK_CLOSE_FRACTION(*root, 0.23, 1e-6);

  BOOST_CHECK_EQUAL(accepted_roots.size(), 0);
  BOOST_CHECK_EQUAL(rejected_roots.size(), 1);
  BOOST_CHECK(!rejected_roots.front().accepted);
  BOOST_CHECK_LT(rejected_roots.front().Z, 0.0);

  BOOST_CHECK_EQUAL(diag.first_interval_shell, 5);
  BOOST_CHECK_EQUAL(diag.first_accepted_shell, -1);
  BOOST_CHECK_EQUAL(diag.chosen_shell, 5);
  BOOST_CHECK_EQUAL(diag.intervals_found, 1);
}

BOOST_AUTO_TEST_SUITE_END()