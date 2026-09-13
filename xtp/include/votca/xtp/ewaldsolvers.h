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
#ifndef VOTCA_XTP_EWALDSOLVERS_H
#define VOTCA_XTP_EWALDSOLVERS_H

// Standard includes
#include <chrono>
#include <cmath>
#include <limits>

// Third party includes
#include <Eigen/Eigenvalues>
#include <boost/format.hpp>

// Local VOTCA includes
#include "votca/tools/constants.h"
#include "votca/xtp/eigen.h"
#include "votca/xtp/ewaldperiodicdipoleoperator.h"
#include "votca/xtp/ewaldregistry.h"
#include "votca/xtp/logger.h"

/**
 * \brief Iterative solvers for the periodic induced-dipole equation
 *        A*mu = F_perm, with A = EwaldPeriodicDipoleOperator.
 *
 * Extracted from EwaldBackground (where they began life as file-local
 * helpers in a calculator header) for two reasons: the MM/MM and QM/MM
 * calculators need the same solvers without including a calculator, and
 * a private calculator header cannot be unit-tested -- the test target
 * does not put xtp/src/libxtp on its include path, so the solvers went
 * untested while every other Ewald layer had its own test.
 *
 * Two solvers are provided, deliberately:
 *
 *   SolveWithJOR   -- Jacobi over-relaxation, reproducing legacy
 *                     PolarBackground's own SOR scheme (same criterion,
 *                     same omega) so the two codes can be compared
 *                     iteration by iteration.
 *   SolveWithIndefinitenessCheck -- preconditioned conjugate gradient,
 *                     far fewer iterations on the same system, with a
 *                     per-iteration curvature check (see below).
 *
 * Both converge to the same fixed point; that equivalence is asserted
 * directly in test_ewaldsolvers.cc rather than assumed.
 */

namespace votca {
namespace xtp {

// Runs the standard preconditioned CG algorithm directly (rather than via
// Eigen::ConjugateGradient) so every iteration's own curvature,
// p.dot(A*p), can be inspected. For a genuinely SPD operator this is
// mathematically guaranteed positive for every nonzero p, at every
// iteration -- CG's own convergence theory depends on it. If it is ever
// <= 0, that is a direct proof (not an inference from residual behavior)
// that the operator is NOT positive-definite, and this stops immediately
// with a clear diagnostic rather than continue computing meaningless
// iterations. This exists specifically to distinguish, for a real
// divergence found this session (residual growing past its own starting
// value with a plain PCG solve, on a system legacy's own SOR handles
// fine at the same physics), between two genuinely different
// explanations that could not be told apart from residual behavior
// alone: the true global operator being indefinite (this check would
// fire), versus a merely poorly-preconditioned but still SPD system
// (this check would never fire, even on a slowly-converging or stalled
// run). Templated on the preconditioner type for the same reason
// EwaldBackground::Evaluate's own PCG solve needs two separate branches
// below: the preconditioner type is a compile-time choice here too, and
// this avoids writing the same loop out twice.
struct PcgIndefinitenessResult {
  Eigen::VectorXd x;
  Index iterations = 0;
  double residual = 0.0;
  bool converged = false;
  // -1 if the operator was never found indefinite; otherwise the
  // 1-based iteration number at which p.dot(A*p) <= 0 was first found.
  // This is a real, unambiguous certificate whenever it fires (CG's own
  // convergence theory strictly requires p.dot(A*p) > 0 for every
  // nonzero p if the operator is genuinely SPD), but it can miss real
  // indefiniteness for a "lucky" (or, from this test's own point of
  // view, unlucky) b that happens not to excite the negative-eigenvalue
  // direction at all -- confirmed directly with a deliberately
  // constructed 3x3 counterexample (a genuinely indefinite matrix,
  // eigenvalues -1/3/5, where b=[1,1,1] is EXACTLY orthogonal to the
  // eigenvector for -1, letting CG converge to the exact right answer
  // in 2 iterations without curvature ever going non-positive). See
  // lanczos_min_eigenvalue below for a real (if still b-dependent, just
  // less narrowly so) improvement on this.
  Index indefinite_at_iteration = -1;
  double indefinite_curvature = 0.0;
  // The smallest eigenvalue of the (small, m x m, m = iterations run)
  // Lanczos tridiagonal matrix built from every alpha_k/beta_k this CG
  // run itself computed -- a standard result (the CG-Lanczos
  // correspondence; verified directly here against known eigenvalues
  // for both an SPD test matrix, exact match to 1e-14, and a genuinely
  // indefinite one with a GENERIC right-hand side, exact match there
  // too) that these Ritz values approximate the extreme eigenvalues of
  // the (preconditioned) operator using the WHOLE accumulated Krylov
  // subspace, not one iteration's own snapshot the way
  // indefinite_at_iteration is. Still not a b-independent guarantee (the
  // same 3x3 counterexample above, with its own adversarially exact
  // orthogonality, defeats this too -- confirmed directly) but
  // meaningfully more robust for any b that is not exactly orthogonal to
  // the relevant eigenvector, which a genuine physical field vector,
  // not an adversarially constructed one, essentially never is exactly.
  // NaN if fewer than 2 iterations ran (need at least 2 alphas for a
  // meaningful 2x2 tridiagonal matrix).
  double lanczos_min_eigenvalue = std::numeric_limits<double>::quiet_NaN();
};

template <typename Preconditioner>
PcgIndefinitenessResult SolveWithIndefinitenessCheck(
    const EwaldPeriodicDipoleOperator& op, Preconditioner& precond,
    const Eigen::VectorXd& b, Index max_iter, double tol, Logger& log,
    std::chrono::steady_clock::time_point t_start) {
  auto elapsed_s = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         t_start)
        .count();
  };

  PcgIndefinitenessResult result;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
  Eigen::VectorXd residual_vec = b - op * x;
  const double rhs_norm2 = b.squaredNorm();
  const double threshold =
      std::max(tol * tol * rhs_norm2, std::numeric_limits<double>::min());
  double residual_norm2 = residual_vec.squaredNorm();
  double tol_error = std::sqrt(residual_norm2 / rhs_norm2);
  bool converged = (residual_norm2 < threshold);
  // Every alpha_k/beta_k this run computes, in order -- see
  // lanczos_min_eigenvalue's own documentation for what these build.
  std::vector<double> alphas;
  std::vector<double> betas;

  if (!converged) {
    Eigen::VectorXd p = precond.solve(residual_vec);
    Eigen::VectorXd z(b.size()), tmp(b.size());
    double abs_new = residual_vec.dot(p);
    Index i = 0;
    while (i < max_iter) {
      tmp.noalias() = op * p;
      const double curvature = p.dot(tmp);
      XTP_LOG(Log::info, log)
          << TimeStamp() << "   PCG iter " << (i + 1)
          << ": curvature p.A.p=" << curvature << " (" << elapsed_s() << "s)"
          << std::flush;
      if (curvature <= 0.0) {
        result.indefinite_at_iteration = i + 1;
        result.indefinite_curvature = curvature;
        XTP_LOG(Log::info, log)
            << TimeStamp() << "   PCG iter " << (i + 1)
            << ": p.A.p = " << curvature
            << " <= 0 -- the operator is NOT positive-definite (this is a "
               "direct algebraic certificate, not an inference from "
               "residual behavior). Stopping here rather than continue "
               "computing iterations that CG's own convergence theory no "
               "longer covers."
            << std::flush;
        break;
      }
      const double alpha = abs_new / curvature;
      alphas.push_back(alpha);
      x += alpha * p;
      residual_vec -= alpha * tmp;

      residual_norm2 = residual_vec.squaredNorm();
      tol_error = std::sqrt(residual_norm2 / rhs_norm2);
      if (residual_norm2 < threshold) {
        converged = true;
        ++i;
        break;
      }

      z = precond.solve(residual_vec);
      const double abs_old = abs_new;
      abs_new = residual_vec.dot(z);
      const double beta = abs_new / abs_old;
      betas.push_back(beta);
      p = z + beta * p;
      ++i;
    }
    result.iterations = i;
  } else {
    result.iterations = 0;
  }

  // Lanczos tridiagonal matrix from alphas/betas -- see
  // lanczos_min_eigenvalue's own documentation for the derivation and
  // its own direct numerical verification (both SPD and indefinite test
  // cases, exact match to the true eigenvalues in both). Needs at least
  // 2 alphas (i.e. 1 beta) for a meaningful 2x2 matrix; a 1x1 "matrix"
  // is just 1/alphas[0], not a real eigenvalue ESTIMATE of anything.
  if (alphas.size() >= 2) {
    const std::size_t m = alphas.size();
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(m, m);
    T(0, 0) = 1.0 / alphas[0];
    for (std::size_t i = 1; i < m; ++i) {
      T(i, i) = 1.0 / alphas[i] + betas[i - 1] / alphas[i - 1];
      const double off = std::sqrt(betas[i - 1]) / alphas[i - 1];
      T(i, i - 1) = off;
      T(i - 1, i) = off;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
    result.lanczos_min_eigenvalue = es.eigenvalues()(0);
  }

  result.x = x;
  result.residual = tol_error;
  result.converged = converged;
  return result;
}

// Per-site polarizability blocks (P, not P^-1), used for the plain
// Jacobi-Over-Relaxation (JOR) iteration below. Confirmed directly
// against legacy PolarBackground's own source (not re-derived from
// theory alone) that this, not EwaldBlockJacobiPreconditioner's own
// larger (intramolecular-Thole-coupled) block, is the exact match to
// legacy's own per-site update: APolarSite::Induce computes
//   mu_new = (1-wSOR)*mu_old + wSOR*(-P)*(F_perm+F_induced)
// which is a textbook weighted-Jacobi step once D = P^-1 is recognized
// as this system's own diagonal block (EwaldPeriodicDipoleOperator's
// own diagonal is built from getPInv() exactly) -- D^-1 = P follows
// directly, with no re-derivation needed of the specific field-based
// expression legacy itself uses internally.
//
// Also confirmed directly (by tracing legacy's own PolarBackground::
// Evaluate) that legacy's own "SOR" is Jacobi-style, not sequential
// Gauss-Seidel-based SOR despite the shared name: legacy computes EVERY
// site's own field in full (step III.B, "(Re-)generate induction
// fields") before updating ANY site's own dipole (step III.C, "Induce
// again") -- every field used in a given iteration is built entirely
// from the PREVIOUS iteration's own induced dipoles.
class EwaldSitePolarizabilityBlocks {
 public:
  EwaldSitePolarizabilityBlocks(const EwaldRegistry& registry,
                                const std::vector<Index>& ids) {
    for (Index id : ids) {
      const PolarSegment& segment =
          registry.Get(id, EwaldChargeState::Neutral);
      for (const PolarSite& site : segment) {
        // getPInv() is this codebase's own stored quantity (see
        // EwaldPeriodicDipoleOperator's own diagonal); P itself is not
        // separately exposed anywhere, so recovering it needs a real
        // 3x3 inverse here -- cheap and exact for a non-singular
        // polarizability tensor.
        blocks_.push_back(site.getPInv().inverse());
      }
    }
  }

  Index size() const { return 3 * Index(blocks_.size()); }

  // Exposes the per-site P block directly.
  const Eigen::Matrix3d& GetBlock(Index n) const { return blocks_[n]; }

  // Returns P * v, applied per 3-DOF site block -- the D^-1 * v the
  // JOR iteration below needs (see this class's own documentation for
  // why D^-1 = P here specifically). Site ordering matches
  // EwaldPeriodicDipoleOperator's own targets_ layout exactly (same
  // registry, same ids, same per-segment site iteration order), so no
  // separate offset bookkeeping is needed here: site n's own 3 DOF are
  // always at [3n, 3n+3) for both this class and that operator.
  Eigen::VectorXd Apply(const Eigen::VectorXd& v) const {
    Eigen::VectorXd result(size());
    for (std::size_t n = 0; n < blocks_.size(); ++n) {
      result.segment<3>(3 * Index(n)) =
          blocks_[n] * v.segment<3>(3 * Index(n));
    }
    return result;
  }

 private:
  std::vector<Eigen::Matrix3d> blocks_;
};

struct JorResult {
  Eigen::VectorXd x;
  Index iterations = 0;
  double residual = 0.0;
  double max_dU = 0.0;
  double avg_dU = 0.0;
  bool converged = false;
};

// Mirrors legacy APolarSite::HistdU() exactly (confirmed directly
// against its own source, apolarsite.cc): per-site relative change in
// induced dipole moment between successive iterates, with a small-
// magnitude fallback so a near-zero dipole (e.g. this site's own first
// iteration, starting from x=0) doesn't divide by ~zero. Legacy's own
// "small" constant is documented there as "1e-10 e*nm" -- this
// codebase uses atomic units (bohr) throughout, not nm (see e.g.
// EwaldRealSpaceSum's own class documentation), so using legacy's raw
// 1e-10 unconverted would apply the wrong ABSOLUTE threshold (1e-10 nm
// and 1e-10 bohr are very different physical dipole magnitudes) --
// tools::conv::nm2bohr is this codebase's own existing conversion
// constant (already used elsewhere in this file for r_min), reused
// here rather than a separately hardcoded factor. epstol itself (see
// SolveWithJOR below) needs no such conversion: it's a ratio of two
// same-unit quantities, genuinely dimensionless regardless of which
// unit system U0/U1 happen to be expressed in.
double SiteRelativeDipoleChange(const Eigen::Vector3d& U0,
                                const Eigen::Vector3d& U1) {
  const double small = 1e-10 * tools::conv::nm2bohr;
  const double abs_U0 = U0.norm();
  const double abs_U1 = U1.norm();
  const double abs_dU = (U1 - U0).norm();
  const double abs_U = (abs_U1 > abs_U0 && small > abs_U0) ? abs_U1 : abs_U0;
  if (small > abs_U) {
    return abs_U;
  }
  return abs_dU / abs_U;
}

// Plain weighted-Jacobi iteration for A*x=b (A = EwaldPeriodicDipoleOperator),
// matching legacy PolarBackground's own SOR exactly in structure (see
// EwaldSitePolarizabilityBlocks's own documentation for the direct
// confirmation against legacy's own source). Iteration, in the general
// weighted-Jacobi form for A = D + R:
//   x_new = x_old + omega * D^-1 * (b - A*x_old)
// Unlike PCG, this makes no positive-definiteness assumption anywhere
// in its own derivation -- its convergence instead depends on the
// spectral radius of the iteration matrix (I - omega*D^-1*A) staying
// below 1, a genuinely different (and empirically, for this kind of
// system, more forgiving) condition. Verified directly, before writing
// any of this, against a small SPD test system (exact match to the
// direct solve) and against a genuinely indefinite one (the same 3x3
// counterexample used to validate this session's own indefiniteness
// checks) -- JOR converges on the indefinite case where CG structurally
// cannot, matching the entire reason for using it here.
//
// Convergence criterion matches legacy PolarBackground's own exactly
// (confirmed directly against its own source, step III.D) -- NOT the
// global relative residual ||b-Ax||/||b|| this function's own first
// version used (borrowed, incorrectly, from the PCG code path's own
// stopping rule): legacy checks a per-site relative change in induced
// dipole moment between successive iterates (SiteRelativeDipoleChange
// above) against epstol=1e-3 (legacy's own hardcoded value, not this
// codebase's own pcg_tolerance_) for EVERY site, with a looser
// avgdU < epstol*0.1 fallback that can force convergence even if a few
// individual sites remain just above epstol. The global residual is
// still computed and logged every iteration (informational only -- it
// is a genuinely different quantity, not proportional to max_dU/avg_dU
// in general, and no longer decides when this function itself stops).
//
// match_legacy_first_step, when true, uses omega=1.0 (fully unrelaxed)
// for iteration 0 only, then the configured omega for every iteration
// after -- literally replicating legacy's own two-phase structure
// (InduceDirect(), fully unrelaxed, THEN the main loop's own wSOR-
// relaxed Induce() calls) rather than applying the same omega-relaxed
// step from x=0 that legacy never does. Added after this session's own
// attempt to algebraically relate the two codebases' differing
// recursions (legacy: mu2 = mu1 - wSOR*P*FU(mu1), coefficient 1 on
// mu1, since legacy's own mu1 is already fully unrelaxed; this
// function's own un-adjusted recursion: x2 = (2-omega)*x1 +
// omega*P*C*x1, a different coefficient, purely because THIS x1 was
// only omega-relaxed to begin with) turned out to rest on an
// unverified assumption (that C, this class's own leftover-from-D
// coupling operator, equals legacy's own FU field) that was never
// independently confirmed -- this flag sidesteps that assumption
// entirely by making the two codebases' own x1/mu1 and x2/mu2 directly,
// literally comparable (same relaxation history on both sides), so a
// comparison needs only the same simple bohr<->nm dipole conversion
// already validated for x1 vs mu1, not any field-specific conversion
// or unverified operator decomposition.
JorResult SolveWithJOR(const EwaldPeriodicDipoleOperator& op,
                       const EwaldSitePolarizabilityBlocks& site_p,
                       const Eigen::VectorXd& b, Index max_iter,
                       double omega, Logger& log,
                       std::chrono::steady_clock::time_point t_start,
                       const EwaldRegistry& registry,
                       const std::vector<Index>& ids,
                       bool match_legacy_first_step) {
  auto elapsed_s = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         t_start)
        .count();
  };

  const double kEpsTol = 1e-3;  // legacy's own hardcoded value

  JorResult result;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
  const double rhs_norm = b.norm();
  const Index n_sites = b.size() / 3;

  Index i = 0;
  for (; i < max_iter; ++i) {
    Eigen::VectorXd residual_vec = b - op * x;
    const double residual_norm = residual_vec.norm();
    const Eigen::VectorXd x_old = x;
    const double this_step_omega =
        (i == 0 && match_legacy_first_step) ? 1.0 : omega;
    x = x_old + this_step_omega * site_p.Apply(residual_vec);



    double max_dU = -1.0;
    double avg_dU = 0.0;
    for (Index n = 0; n < n_sites; ++n) {
      const double dU = SiteRelativeDipoleChange(x_old.segment<3>(3 * n),
                                                 x.segment<3>(3 * n));
      avg_dU += dU;
      if (dU > max_dU) {
        max_dU = dU;
      }
    }
    avg_dU /= double(n_sites);

    XTP_LOG(Log::info, log)
        << TimeStamp() << "   JOR iter " << (i + 1) << ": max_dU=" << max_dU
        << " avg_dU=" << avg_dU << " (residual=" << residual_norm / rhs_norm
        << ", informational only -- see this function's own documentation "
           "for why max_dU/avg_dU, not this, is the actual stopping "
           "criterion) ("
        << elapsed_s() << "s)" << std::flush;

    result.residual = residual_norm / rhs_norm;
    result.max_dU = max_dU;
    result.avg_dU = avg_dU;

    bool converged = (max_dU <= kEpsTol);
    if (avg_dU < kEpsTol * 0.1) {
      converged = true;
    }
    if (converged) {
      result.converged = true;
      ++i;
      break;
    }
  }
  result.iterations = i;
  result.x = x;
  return result;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDSOLVERS_H
