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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE ewaldsolvers_test

// Standard includes
#include <iostream>

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/ewaldblockjacobipreconditioner.h"
#include "votca/xtp/ewaldperiodicdipoleoperator.h"
#include "votca/xtp/ewaldrealspacesum.h"
#include "votca/xtp/ewaldreciprocalspacesum.h"
#include "votca/xtp/ewaldregistry.h"
#include "votca/xtp/ewaldshapecorrection.h"
#include "votca/xtp/ewaldsolvers.h"
#include "votca/xtp/logger.h"

// NOTE ON WHAT THIS FILE DELIBERATELY DOES NOT DO.
//
// None of the checks below compare against a hand-assembled reference
// matrix or a stored expected vector. That is deliberate, and it is the
// direct lesson of a real bug: the operator was built as P^-1 + C
// instead of P^-1 - C -- inverting the induced-induced feedback and so
// converging to the wrong dipoles -- and the existing operator tests did
// NOT catch it, because their dense references had been written by
// mirroring the implementation's own sign convention. They agreed with
// the bug. The same happened a second time with the self-field term.
//
// So every assertion here is either an invariant the true solution must
// satisfy regardless of how it is computed (A*x = b), an agreement
// between two independent solvers, or a structural property of the
// operator (symmetry, positive-definiteness, conditioning). None of
// those can be satisfied by a reference that mirrors a mistake.

using namespace votca;
using namespace votca::xtp;

namespace {

// A small periodic system with several multi-site segments. Deliberately
// multi-site: the intramolecular coupling block is where two of the
// historical sign mistakes lived, and a single-site system has none.
EwaldRegistry MakeSystem(double box_length, Index n_per_side) {
  EwaldRegistry registry;
  const double d = box_length / double(n_per_side);
  Index id = 0;
  for (Index a = 0; a < n_per_side; ++a) {
    for (Index b = 0; b < n_per_side; ++b) {
      for (Index c = 0; c < n_per_side; ++c) {
        const Eigen::Vector3d centre(double(a) * d, double(b) * d,
                                     double(c) * d);
        PolarSegment seg("seg", id);
        // Tetrahedral, methane-like: four satellites around a centre.
        const double t = 0.63;
        const Eigen::Vector3d offsets[5] = {
            {0.0, 0.0, 0.0}, {t, t, t}, {t, -t, -t}, {-t, t, -t}, {-t, -t, t}};
        for (Index s = 0; s < 5; ++s) {
          PolarSite site(s, (s == 0) ? "C" : "H", centre + offsets[s]);
          site.setpolarization(((s == 0) ? 8.0 : 3.0) *
                               Eigen::Matrix3d::Identity());
          site.setCharge((s == 0) ? -0.4 : 0.1);
          site.setStaticDipole(
              Eigen::Vector3d(0.01 * double(s + 1), -0.005, 0.002 * double(s)));
          seg.push_back(site);
        }
        registry.Register(id, EwaldChargeState::Neutral, seg);
        ++id;
      }
    }
  }
  return registry;
}

std::vector<Index> AllIds(Index n) {
  std::vector<Index> ids;
  for (Index i = 0; i < n; ++i) {
    ids.push_back(i);
  }
  return ids;
}

// Right-hand side: any fixed, reproducible vector will do. The solvers
// are being tested, not the permanent-field generation (which has its
// own tests), so this deliberately does not depend on that machinery.
Eigen::VectorXd MakeRhs(Index size) {
  Eigen::VectorXd b(size);
  for (Index i = 0; i < size; ++i) {
    b(i) = 1e-3 * std::sin(0.7 * double(i) + 0.3);
  }
  return b;
}

// Dense form of the operator, obtained by applying it to unit vectors.
// This is a MEASUREMENT of the operator, not a re-derivation of it, so
// it cannot encode a mistaken sign convention of its own.
Eigen::MatrixXd DenseOperator(const EwaldPeriodicDipoleOperator& op) {
  const Index n = op.rows();
  Eigen::MatrixXd A(n, n);
  for (Index i = 0; i < n; ++i) {
    Eigen::VectorXd e = Eigen::VectorXd::Zero(n);
    e(i) = 1.0;
    A.col(i) = op * e;
  }
  return A;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(ewaldsolvers_test)

// The operator must be symmetric positive definite: CG's convergence
// theory requires it, and the block-Jacobi preconditioner's own LDLT
// factorization assumes it. This is also the cheapest possible guard on
// the coupling sign -- flipping C from -C to +C breaks definiteness long
// before it breaks anything a casual eyeball check would notice.
BOOST_AUTO_TEST_CASE(operator_is_symmetric_positive_definite) {
  const double alpha = 0.3;
  const double L = 14.0;
  Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  EwaldRegistry registry = MakeSystem(L, 2);
  std::vector<Index> ids = AllIds(8);

  // Tighter real-space convergence than the other cases here, on
  // purpose. EwaldRealSpaceSum converges each target's own neighbour
  // list independently, so at looser settings a pair can be summed for
  // target i and not for target j, which makes A asymmetric at the level
  // of the truncation itself -- measured at 8e-5 with r_min=12 /
  // field_tol=1e-12, dropping to 1e-16 with the values below. That is a
  // property of the adaptive shell search, not a bug, but it means
  // symmetry can only be asserted at machine precision once the sum is
  // converged tightly enough that the neighbour sets agree.
  EwaldRealSpaceSum real_sum(box, registry, alpha, 0.39, /*r_min=*/30.0,
                             /*field_tol=*/1e-14);
  EwaldReciprocalSpaceSum recip_sum(box, registry, alpha, /*k_max=*/8.0);
  EwaldShapeCorrection shape(box.determinant(), registry, EwaldShape::Cube);
  EwaldPeriodicDipoleOperator op(registry, real_sum, recip_sum, shape, ids,
                                 alpha, 0.39);

  const Eigen::MatrixXd A = DenseOperator(op);

  const double asymmetry = (A - A.transpose()).norm() / A.norm();
  BOOST_CHECK_SMALL(asymmetry, 1e-12);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(0.5 * (A + A.transpose()));
  BOOST_REQUIRE(es.info() == Eigen::Success);
  const double lambda_min = es.eigenvalues().minCoeff();
  std::cout << "asymmetry = " << asymmetry << "  lambda_min = " << lambda_min
            << "  lambda_max = " << es.eigenvalues().maxCoeff() << std::endl;
  BOOST_CHECK_GT(lambda_min, 0.0);
}

// The defining property of the answer, independent of how it was found:
// the returned x must actually solve A*x = b. A stored reference vector
// would only assert that the code still does what it did before, which
// is exactly what let an inverted coupling sign survive.
BOOST_AUTO_TEST_CASE(pcg_solution_satisfies_the_system) {
  const double alpha = 0.3;
  const double L = 14.0;
  Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  EwaldRegistry registry = MakeSystem(L, 2);
  std::vector<Index> ids = AllIds(8);

  EwaldRealSpaceSum real_sum(box, registry, alpha, 0.39, 12.0, 1e-12);
  EwaldReciprocalSpaceSum recip_sum(box, registry, alpha, 8.0);
  EwaldShapeCorrection shape(box.determinant(), registry, EwaldShape::Cube);
  EwaldPeriodicDipoleOperator op(registry, real_sum, recip_sum, shape, ids,
                                 alpha, 0.39);

  Logger log;
  log.setReportLevel(Log::error);
  const Eigen::VectorXd b = MakeRhs(op.rows());

  EwaldBlockJacobiPreconditioner precond(registry, ids, alpha, 0.39);
  auto result = SolveWithIndefinitenessCheck(op, precond, b, /*max_iter=*/200,
                                             /*tol=*/1e-12, log,
                                             std::chrono::steady_clock::now());

  BOOST_REQUIRE(result.converged);
  // Never found indefinite: a direct certificate, not an inference.
  BOOST_CHECK_EQUAL(result.indefinite_at_iteration, -1);

  const double rel_residual = (op * result.x - b).norm() / b.norm();
  std::cout << "PCG iterations = " << result.iterations
            << "  relative residual = " << rel_residual << std::endl;
  BOOST_CHECK_SMALL(rel_residual, 1e-10);
  BOOST_CHECK_GT(result.x.norm(), 1e-10);  // not vacuously zero
}

// JOR and PCG walk entirely different paths through the same space. If
// both are correct they must land on the same fixed point; if either has
// a solver-specific mistake they will not. Neither side of this
// comparison is a hand-written reference, so neither can encode the
// other's error.
BOOST_AUTO_TEST_CASE(jor_and_pcg_agree) {
  const double alpha = 0.3;
  const double L = 14.0;
  Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  EwaldRegistry registry = MakeSystem(L, 2);
  std::vector<Index> ids = AllIds(8);

  EwaldRealSpaceSum real_sum(box, registry, alpha, 0.39, 12.0, 1e-12);
  EwaldReciprocalSpaceSum recip_sum(box, registry, alpha, 8.0);
  EwaldShapeCorrection shape(box.determinant(), registry, EwaldShape::Cube);
  EwaldPeriodicDipoleOperator op(registry, real_sum, recip_sum, shape, ids,
                                 alpha, 0.39);

  Logger log;
  log.setReportLevel(Log::error);
  const Eigen::VectorXd b = MakeRhs(op.rows());
  const auto t0 = std::chrono::steady_clock::now();

  EwaldBlockJacobiPreconditioner precond(registry, ids, alpha, 0.39);
  auto pcg = SolveWithIndefinitenessCheck(op, precond, b, 500, 1e-12, log, t0);
  BOOST_REQUIRE(pcg.converged);

  EwaldSitePolarizabilityBlocks site_p(registry, ids);
  auto jor = SolveWithJOR(op, site_p, b, /*max_iter=*/20000, /*omega=*/0.35,
                          log, t0, /*match_legacy_first_step=*/false);

  // JOR stops on a per-site dipole-change criterion with a hardcoded
  // 1e-3 threshold (matching legacy's own epstol), NOT on a residual, so
  // it cannot be driven to the precision PCG reaches here however many
  // iterations it is given. Demanding agreement at PCG's tolerance would
  // therefore fail for a reason that has nothing to do with correctness.
  //
  // The meaningful invariant is that the gap between the two is
  // explained by JOR's own incomplete convergence: if the solvers had
  // genuinely different fixed points, the difference between them would
  // NOT be of the order of JOR's own residual, it would be some
  // unrelated constant.
  const double jor_residual = (op * jor.x - b).norm() / b.norm();
  const double rel_diff = (jor.x - pcg.x).norm() / pcg.x.norm();
  std::cout << "JOR iterations = " << jor.iterations
            << "  JOR residual = " << jor_residual
            << "  |x_JOR - x_PCG|/|x_PCG| = " << rel_diff << std::endl;

  BOOST_REQUIRE(jor.converged);
  BOOST_REQUIRE_SMALL(jor_residual, 1e-2);
  BOOST_CHECK_LT(rel_diff, 10.0 * jor_residual);
}

// A preconditioner cannot change the fixed point, only the path to it,
// so no correctness test can detect a broken one -- it just silently
// costs iterations. This is the guard for that: a real bug once left the
// block-Jacobi coupling term with the wrong sign, which made the
// conditioning ~4x WORSE than no preconditioner at all and nearly
// doubled the iteration count, while every result stayed correct.
BOOST_AUTO_TEST_CASE(preconditioner_improves_conditioning) {
  const double alpha = 0.3;
  const double L = 14.0;
  Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  EwaldRegistry registry = MakeSystem(L, 2);
  std::vector<Index> ids = AllIds(8);

  EwaldRealSpaceSum real_sum(box, registry, alpha, 0.39, 12.0, 1e-12);
  EwaldReciprocalSpaceSum recip_sum(box, registry, alpha, 8.0);
  EwaldShapeCorrection shape(box.determinant(), registry, EwaldShape::Cube);
  EwaldPeriodicDipoleOperator op(registry, real_sum, recip_sum, shape, ids,
                                 alpha, 0.39);

  const Eigen::MatrixXd A = DenseOperator(op);
  const Index n = A.rows();

  // Dense form of the preconditioner, again by measurement.
  EwaldBlockJacobiPreconditioner precond(registry, ids, alpha, 0.39);
  Eigen::MatrixXd M(n, n);
  for (Index i = 0; i < n; ++i) {
    Eigen::VectorXd e = Eigen::VectorXd::Zero(n);
    e(i) = 1.0;
    M.col(i) = precond.solve(e);
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esA(0.5 * (A + A.transpose()));
  BOOST_REQUIRE(esA.info() == Eigen::Success);
  const double cond_A =
      esA.eigenvalues().maxCoeff() / esA.eigenvalues().minCoeff();

  const Eigen::MatrixXd MA = M * A;
  Eigen::EigenSolver<Eigen::MatrixXd> esMA(MA);
  BOOST_REQUIRE(esMA.info() == Eigen::Success);
  double lo = std::numeric_limits<double>::max();
  double hi = -std::numeric_limits<double>::max();
  for (Index i = 0; i < n; ++i) {
    const double re = esMA.eigenvalues()(i).real();
    lo = std::min(lo, re);
    hi = std::max(hi, re);
  }
  BOOST_REQUIRE_GT(lo, 0.0);
  const double cond_MA = hi / lo;

  std::cout << "cond(A) = " << cond_A << "  cond(M*A) = " << cond_MA
            << std::endl;
  BOOST_CHECK_LT(cond_MA, cond_A);
}

// The operator parallelizes over targets with no cross-target reduction,
// so the summation order is fixed and the result must be independent of
// thread count -- not merely close, but bit-for-bit identical. Asserting
// exact equality is deliberate: a tolerance here would quietly permit a
// genuine race to creep in.
BOOST_AUTO_TEST_CASE(operator_is_thread_count_independent) {
  const double alpha = 0.3;
  const double L = 14.0;
  Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  EwaldRegistry registry = MakeSystem(L, 2);
  std::vector<Index> ids = AllIds(8);

  EwaldRealSpaceSum real_sum(box, registry, alpha, 0.39, 12.0, 1e-12);
  EwaldReciprocalSpaceSum recip_sum(box, registry, alpha, 8.0);
  EwaldShapeCorrection shape(box.determinant(), registry, EwaldShape::Cube);
  EwaldPeriodicDipoleOperator op(registry, real_sum, recip_sum, shape, ids,
                                 alpha, 0.39);

  const Eigen::VectorXd v = MakeRhs(op.rows());

#ifdef _OPENMP
  const Index saved_threads = Index(omp_get_max_threads());
  omp_set_num_threads(1);
#endif
  const Eigen::VectorXd serial = op * v;
#ifdef _OPENMP
  omp_set_num_threads(4);
#endif
  const Eigen::VectorXd parallel = op * v;
#ifdef _OPENMP
  omp_set_num_threads(int(saved_threads));
#endif

  BOOST_REQUIRE_EQUAL(serial.size(), parallel.size());
  bool identical = true;
  for (Index i = 0; i < serial.size(); ++i) {
    if (serial(i) != parallel(i)) {
      identical = false;
      break;
    }
  }
  BOOST_CHECK(identical);
}

// The real-space distance cull discards pairs whose erfc-screened
// contribution is below the working precision. It is a pure cost
// optimization: at the default screening_factor the field must be
// unchanged. This guards it as alpha varies -- the cutoff is
// screening_factor/alpha, so a system with different alpha gets a
// different cutoff without anyone touching the option.
BOOST_AUTO_TEST_CASE(screening_cutoff_does_not_change_the_field) {
  const double alpha = 0.3;
  const double L = 14.0;
  Eigen::Matrix3d box = L * Eigen::Matrix3d::Identity();
  EwaldRegistry registry = MakeSystem(L, 2);

  // Default cutoff versus an effectively infinite one.
  EwaldRealSpaceSum culled(box, registry, alpha, 0.39, 12.0, 1e-12,
                           /*shell_width=*/0.945, /*n_max=*/15,
                           /*screening_factor=*/6.0);
  EwaldRealSpaceSum unculled(box, registry, alpha, 0.39, 12.0, 1e-12,
                             /*shell_width=*/0.945, /*n_max=*/15,
                             /*screening_factor=*/1e9);

  double worst = 0.0;
  for (Index id = 0; id < 8; ++id) {
    PolarSite a = registry.Get(id, EwaldChargeState::Neutral)[0];
    PolarSite b = a;
    a.Reset();
    b.Reset();
    culled.AddFieldAt<Estatic::V>(id, a, EwaldChargeState::Neutral);
    unculled.AddFieldAt<Estatic::V>(id, b, EwaldChargeState::Neutral);
    const double denom = std::max(b.V().norm(), 1e-30);
    worst = std::max(worst, (a.V() - b.V()).norm() / denom);
  }
  std::cout << "worst relative field change from the cull = " << worst
            << std::endl;
  BOOST_CHECK_SMALL(worst, 1e-12);

  const auto stats = culled.GetNeighborStats();
  std::cout << "culled fraction = " << stats.culled_fraction() << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
