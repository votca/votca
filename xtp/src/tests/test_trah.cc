/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

#define BOOST_TEST_MODULE trah_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/trah.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>

using namespace votca::xtp;
using namespace votca;

namespace {
// Deterministic pseudo-random numbers in [-1, 1) from a 64-bit linear
// congruential generator: the same with every Eigen version and C library
// (Eigen::Random is not). Seeded per test via Seed(); all tests were
// checked to pass with the seeds offset by 0-40.
std::uint64_t random_state = 1;
void Seed(std::uint64_t seed) { random_state = seed; }
double NextRandom() {
  random_state = random_state * 6364136223846793005ULL + 1442695040888963407ULL;
  return double(random_state >> 11) * (2.0 / 9007199254740992.0) - 1.0;
}
Eigen::MatrixXd RandomMatrix(Index rows, Index cols) {
  Eigen::MatrixXd m(rows, cols);
  for (Index j = 0; j < cols; ++j) {
    for (Index i = 0; i < rows; ++i) {
      m(i, j) = NextRandom();
    }
  }
  return m;
}
Eigen::VectorXd RandomVector(Index n) { return RandomMatrix(n, 1).col(0); }

Eigen::MatrixXd RandomSymmetric(Index n, double shift) {
  Eigen::MatrixXd A = RandomMatrix(n, n);
  return 0.5 * (A + A.transpose()) + shift * Eigen::MatrixXd::Identity(n, n);
}

trah::HessianProduct Dense(const Eigen::MatrixXd& H, Index& calls) {
  return [&H, &calls](const Eigen::MatrixXd& V) {
    calls += V.cols();
    return Eigen::MatrixXd(H * V);
  };
}

double Model(const Eigen::VectorXd& g, const Eigen::MatrixXd& H,
             const Eigen::VectorXd& x) {
  return g.dot(x) + 0.5 * x.dot(H * x);
}

// exp of an antisymmetric matrix by scaling and squaring of a Taylor series
Eigen::MatrixXd ExpAntisymmetric(const Eigen::MatrixXd& K) {
  const Index n = K.rows();
  const int squarings = 10;
  const Eigen::MatrixXd A = K / std::pow(2.0, squarings);
  Eigen::MatrixXd term = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd result = term;
  for (int k = 1; k < 20; ++k) {
    term = (term * A / double(k)).eval();
    result += term;
  }
  for (int k = 0; k < squarings; ++k) {
    result = (result * result).eval();
  }
  return result;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(trah_test)

BOOST_AUTO_TEST_CASE(step_inside_trust_radius_solves_shifted_newton) {
  // Large trust radius: a = 1, and the step solves (H - mu) x = -g with the
  // mu of the dense augmented Hessian.
  Seed(5);
  const Index n = 40;
  const Eigen::MatrixXd H = RandomSymmetric(n, 6.0);
  const Eigen::VectorXd g = 0.1 * RandomVector(n);
  const Eigen::VectorXd diag = H.diagonal();
  Index calls = 0;
  trah::StepOptions opt;
  opt.max_forcing = 1e-8;
  trah::StepResult res = trah::Step(g, diag, Dense(H, calls), 10.0, opt);
  BOOST_CHECK(res.converged);
  BOOST_CHECK_EQUAL(res.alpha, 1.0);

  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n + 1, n + 1);
  M.block(0, 1, 1, n) = g.transpose();
  M.block(1, 0, n, 1) = g;
  M.bottomRightCorner(n, n) = H;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
  const Eigen::VectorXd y = es.eigenvectors().col(0);
  const Eigen::VectorXd x_ref = y.tail(n) / y(0);
  BOOST_CHECK_SMALL((res.x - x_ref).norm() / x_ref.norm(), 1e-6);
  BOOST_CHECK_CLOSE(res.predicted, Model(g, H, res.x), 1e-6);
  BOOST_CHECK_LT(calls, n);
}

BOOST_AUTO_TEST_CASE(step_on_trust_radius) {
  // Small trust radius: |x| within [0.8, 1] of it, and the model decrease
  // at least that of steepest descent of the same length (the Cauchy
  // point, which the first micro-iteration already gives).
  Seed(6);
  const Index n = 60;
  const Eigen::MatrixXd H = RandomSymmetric(n, 5.0);
  const Eigen::VectorXd g = RandomVector(n);
  Index calls = 0;
  const double trust = 0.05;
  trah::StepResult res = trah::Step(g, H.diagonal(), Dense(H, calls), trust);
  BOOST_CHECK(res.converged);
  BOOST_CHECK_GT(res.alpha, 1.0);
  BOOST_CHECK_LE(res.x.norm(), trust * (1 + 1e-12));
  BOOST_CHECK_GE(res.x.norm(), 0.8 * trust);
  const Eigen::VectorXd sd = -res.x.norm() * g.normalized();
  BOOST_CHECK_LE(Model(g, H, res.x), Model(g, H, sd) + 1e-12);
  BOOST_CHECK_CLOSE(res.predicted, Model(g, H, res.x), 1e-6);
}

BOOST_AUTO_TEST_CASE(step_with_negative_curvature) {
  // Indefinite Hessian: the step must still lower the model.
  Seed(7);
  const Index n = 30;
  const Eigen::MatrixXd H = RandomSymmetric(n, 0.0);
  const Eigen::VectorXd g = RandomVector(n);
  Eigen::VectorXd diag = H.diagonal().cwiseAbs().array() + 1.0;
  Index calls = 0;
  trah::StepResult res = trah::Step(g, diag, Dense(H, calls), 0.3);
  BOOST_CHECK_LE(res.x.norm(), 0.3 * (1 + 1e-12));
  BOOST_CHECK_LT(res.predicted, 0.0);
  BOOST_CHECK_CLOSE(res.predicted, Model(g, H, res.x), 1e-6);
}

BOOST_AUTO_TEST_CASE(subspace_restart) {
  // A small subspace limit forces restarts; the step still converges.
  Seed(8);
  const Index n = 80;
  const Eigen::MatrixXd H = RandomSymmetric(n, 8.0);
  const Eigen::VectorXd g = RandomVector(n);
  trah::StepOptions opt;
  opt.max_subspace = 6;
  opt.max_iterations = 200;
  Index calls = 0;
  trah::StepResult res = trah::Step(g, H.diagonal(), Dense(H, calls), 0.2, opt);
  BOOST_CHECK(res.converged);
  BOOST_CHECK_LE(res.x.norm(), 0.2 * (1 + 1e-12));
}

BOOST_AUTO_TEST_CASE(subspace_reuse_after_rejection) {
  // A step repeated with a smaller trust radius for the same gradient
  // starts from the stored subspace: few or no new Hessian products, and
  // the same step as a fresh solve.
  Seed(10);
  const Index n = 70;
  const Eigen::MatrixXd H = RandomSymmetric(n, 3.0);
  const Eigen::VectorXd g = RandomVector(n);
  trah::StepOptions opt;
  opt.max_forcing = 1e-6;  // tight, so that a step needs many products
  trah::Subspace subspace;
  Index calls = 0;
  trah::StepResult first =
      trah::Step(g, H.diagonal(), Dense(H, calls), 0.5, opt, &subspace);
  BOOST_CHECK(!first.reused);
  BOOST_CHECK_EQUAL(subspace.B.cols(), first.products);

  calls = 0;
  trah::StepResult again =
      trah::Step(g, H.diagonal(), Dense(H, calls), 0.1, opt, &subspace);
  BOOST_CHECK(again.reused);
  Index fresh_calls = 0;
  trah::StepResult fresh =
      trah::Step(g, H.diagonal(), Dense(H, fresh_calls), 0.1, opt);
  BOOST_TEST_MESSAGE("first " << first.products << " products, repeated "
                              << again.products << ", fresh "
                              << fresh.products);
  BOOST_CHECK(again.converged);
  BOOST_CHECK_LT(again.products, fresh.products);
  BOOST_CHECK_LE(again.x.norm(), 0.1 * (1 + 1e-12));
  BOOST_CHECK_SMALL(again.predicted - Model(g, H, again.x), 1e-10);
  BOOST_CHECK_SMALL((again.x - fresh.x).norm(), 1e-4 * fresh.x.norm());

  // a different gradient does not reuse
  trah::StepResult other =
      trah::Step(g + 0.1 * Eigen::VectorXd::Ones(n), H.diagonal(),
                 Dense(H, calls), 0.1, opt, &subspace);
  BOOST_CHECK(!other.reused);
}

BOOST_AUTO_TEST_CASE(failed_projected_solve_keeps_last_step) {
  // A direction of strong negative curvature that the gradient almost
  // misses: the first projected problem is solved, the second, after that
  // direction has entered the subspace, has a lowest eigenvector with no
  // component along the augmented coordinate and cannot give a step. The
  // step of the first solve must be returned, consistent with the model
  // (previously the predicted change mixed vectors from the two subspaces).
  Eigen::MatrixXd H(2, 2);
  H << 1.0, 1e-14, 1e-14, -100.0;
  Eigen::VectorXd g(2);
  g << 1.0, 0.0;
  Eigen::VectorXd diag(2);
  diag << 1.0, 100.0;
  trah::StepOptions opt;
  opt.max_forcing = 1e-30;
  opt.residual_floor = 0.0;
  Index calls = 0;
  trah::StepResult res = trah::Step(g, diag, Dense(H, calls), 10.0, opt);
  BOOST_TEST_MESSAGE(
      "products " << res.products << ", iterations " << res.iterations << ", x "
                  << res.x.transpose() << ", predicted " << res.predicted);
  BOOST_CHECK(res.x.allFinite());
  BOOST_CHECK_LT(res.predicted, 0.0);
  BOOST_CHECK_SMALL(res.predicted - Model(g, H, res.x), 1e-12);
  BOOST_CHECK_CLOSE(res.linear, g.dot(res.x), 1e-10);
}

BOOST_AUTO_TEST_CASE(rotate_orbitals_is_exact_exponential) {
  Seed(9);
  const Index n = 9;
  const Index nocc = 4;
  const Index nvirt = n - nocc;
  Eigen::MatrixXd A = RandomMatrix(n, n);
  const Eigen::MatrixXd S =
      A * A.transpose() / double(n) + Eigen::MatrixXd::Identity(n, n);
  // S-orthonormal orbitals
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
  const Eigen::MatrixXd C = es.operatorInverseSqrt();
  const Eigen::VectorXd x = 0.4 * RandomVector(nocc * nvirt);

  Eigen::MatrixXd kappa = Eigen::MatrixXd::Zero(n, n);
  for (Index i = 0; i < nocc; ++i) {
    for (Index a = 0; a < nvirt; ++a) {
      kappa(i, nocc + a) = x(i * nvirt + a);
      kappa(nocc + a, i) = -x(i * nvirt + a);
    }
  }
  const Eigen::MatrixXd ref = C * ExpAntisymmetric(kappa);
  const Eigen::MatrixXd C_new = trah::RotateOrbitals(C, x, nocc);
  BOOST_CHECK_SMALL((C_new - ref).cwiseAbs().maxCoeff(), 1e-10);
  BOOST_CHECK_SMALL(
      (C_new.transpose() * S * C_new - Eigen::MatrixXd::Identity(n, n))
          .cwiseAbs()
          .maxCoeff(),
      1e-10);
  const Eigen::VectorXd ov = trah::OccVirt(kappa, nocc);
  BOOST_CHECK_SMALL((ov - x).norm(), 1e-14);
}

BOOST_AUTO_TEST_CASE(trust_radius_update) {
  trah::TrustRadius t(0.5);
  BOOST_CHECK(t.Update(-1.0, -0.9));  // good model: grow
  BOOST_CHECK_CLOSE(t.Radius(), 0.6, 1e-12);
  BOOST_CHECK(t.Update(-1.0, -0.5));  // keep
  BOOST_CHECK_CLOSE(t.Radius(), 0.6, 1e-12);
  BOOST_CHECK(t.Update(-1.0, -0.1));  // accept, shrink
  BOOST_CHECK_CLOSE(t.Radius(), 0.42, 1e-12);
  BOOST_CHECK(!t.Update(-1.0, 0.1));  // energy rose: reject, shrink
  BOOST_CHECK_CLOSE(t.Radius(), 0.294, 1e-12);
  BOOST_CHECK(t.Update(-1e-12, 2e-12));  // both below the noise floor

  // rejection with the energy along the step: E(t) = linear t + c t^2,
  // c = actual - linear = 0.3, minimum at t* = 0.1 / 0.6 = 1/6
  trah::TrustRadius u(0.8);
  BOOST_CHECK(!u.Update(-0.05, 0.2, -0.1, 0.6));
  BOOST_CHECK_CLOSE(u.Radius(), 0.6 / 6.0, 1e-10);
  // t* outside [0.1, 0.7] is clamped
  trah::TrustRadius v(0.8);
  BOOST_CHECK(!v.Update(-0.05, 10.0, -0.1, 0.5));
  BOOST_CHECK_CLOSE(v.Radius(), 0.05, 1e-10);
  // accepted steps follow the plain update
  trah::TrustRadius w(0.5);
  BOOST_CHECK(w.Update(-1.0, -0.9, -1.5, 0.5));
  BOOST_CHECK_CLOSE(w.Radius(), 0.6, 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
