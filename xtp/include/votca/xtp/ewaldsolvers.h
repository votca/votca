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
 * Kept in a public header (rather than in the EwaldBackground calculator)
 * so that both the MM/MM and QM/MM calculators can use them and so that
 * they can be unit-tested.
 *
 *   SolveWithJOR   -- Jacobi over-relaxation, reproducing legacy
 *                     PolarBackground's SOR scheme (same criterion,
 *                     same omega) for iteration-by-iteration comparison.
 *   SolveWithIndefinitenessCheck -- preconditioned conjugate gradient
 *                     with a per-iteration curvature check.
 *
 * Both converge to the same fixed point (tested in test_ewaldsolvers.cc).
 */

namespace votca {
namespace xtp {

// Preconditioned CG written out explicitly (rather than via
// Eigen::ConjugateGradient) so that the curvature p.dot(A*p) of every
// search direction can be inspected. For an SPD operator it is positive
// for every nonzero p; a value <= 0 proves the operator is not positive
// definite, and the solve stops with a diagnostic. This distinguishes a
// truly indefinite operator from an SPD but poorly preconditioned one,
// which cannot be told apart from the residual history alone.
// Templated on the preconditioner, which is a compile-time choice.
struct PcgIndefinitenessResult {
  Eigen::VectorXd x;
  Index iterations = 0;
  double residual = 0.0;
  bool converged = false;
  // -1 if no non-positive curvature was found; otherwise the 1-based
  // iteration at which p.dot(A*p) <= 0 first occurred. A positive result
  // is a certificate of indefiniteness, but the check can miss it if b
  // has no component along the negative-eigenvalue direction (e.g. the
  // 3x3 matrix with eigenvalues -1/3/5 and b orthogonal to the eigenvector
  // of -1, on which CG converges in 2 iterations).
  Index indefinite_at_iteration = -1;
  double indefinite_curvature = 0.0;
  // Smallest eigenvalue of the m x m Lanczos tridiagonal matrix built from
  // the CG coefficients alpha_k/beta_k (m = iterations run). By the
  // CG-Lanczos correspondence these Ritz values approximate the extreme
  // eigenvalues of the preconditioned operator over the whole Krylov
  // subspace, which is more robust than the per-iteration curvature test,
  // though still b-dependent (it fails on the same orthogonal-b example).
  // NaN if fewer than 2 iterations ran.
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
  // rhs = 0 (e.g. a background with all permanent moments zero) is solved
  // exactly by x = 0: residual_norm2 = 0 is below the threshold DBL_MIN,
  // so the solve returns with zero iterations. The relative residual
  // would be 0/0 and is reported as 0 rather than NaN, which would
  // otherwise make a correct solve look failed in the log. The loop
  // below never runs in that case, so it needs no such guard.
  double tol_error =
      (rhs_norm2 > 0.0) ? std::sqrt(residual_norm2 / rhs_norm2) : 0.0;
  bool converged = (residual_norm2 < threshold);
  // CG coefficients, kept for the Lanczos estimate (lanczos_min_eigenvalue).
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
      XTP_LOG(Log::info, log) << TimeStamp() << "   PCG iter " << (i + 1)
                              << ": curvature p.A.p=" << curvature << " ("
                              << elapsed_s() << "s)" << std::flush;
      if (curvature <= 0.0) {
        result.indefinite_at_iteration = i + 1;
        result.indefinite_curvature = curvature;
        XTP_LOG(Log::info, log)
            << TimeStamp() << "   PCG iter " << (i + 1)
            << ": p.A.p = " << curvature
            << " <= 0 -- the operator is NOT positive-definite. Stopping "
               "CG."
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

  // Lanczos tridiagonal matrix from the CG coefficients:
  //   T_00 = 1/alpha_0, T_ii = 1/alpha_i + beta_{i-1}/alpha_{i-1},
  //   T_{i,i-1} = sqrt(beta_{i-1})/alpha_{i-1}.
  // Requires at least 2 alphas; a 1x1 matrix is not a useful estimate.
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
// Jacobi-Over-Relaxation (JOR) iteration below. These per-site blocks
// (not the larger intramolecular-Thole-coupled blocks of
// EwaldBlockJacobiPreconditioner) reproduce legacy APolarSite::Induce,
//   mu_new = (1-wSOR)*mu_old + wSOR*(-P)*(F_perm+F_induced),
// which is a weighted-Jacobi step with diagonal block D = P^-1 (the
// diagonal of EwaldPeriodicDipoleOperator is built from getPInv()), so
// D^-1 = P.
//
// Despite its name, legacy PolarBackground's "SOR" is Jacobi-style, not
// Gauss-Seidel: all fields are computed from the previous iterate (step
// III.B) before any dipole is updated (step III.C).
class EwaldSitePolarizabilityBlocks {
 public:
  EwaldSitePolarizabilityBlocks(const EwaldRegistry& registry,
                                const std::vector<Index>& ids) {
    for (Index id : ids) {
      const PolarSegment& segment = registry.Get(id, EwaldChargeState::Neutral);
      for (const PolarSite& site : segment) {
        // Only P^-1 is stored on PolarSite; invert the 3x3 block to
        // recover P.
        blocks_.push_back(site.getPInv().inverse());
      }
    }
  }

  Index size() const { return 3 * Index(blocks_.size()); }

  // Exposes the per-site P block directly.
  const Eigen::Matrix3d& GetBlock(Index n) const { return blocks_[n]; }

  // Returns P * v (= D^-1 * v) per 3-DOF site block. Site ordering
  // matches the targets_ layout of EwaldPeriodicDipoleOperator (same
  // registry, ids and site order), so site n occupies [3n, 3n+3) in both.
  Eigen::VectorXd Apply(const Eigen::VectorXd& v) const {
    Eigen::VectorXd result(size());
    for (std::size_t n = 0; n < blocks_.size(); ++n) {
      result.segment<3>(3 * Index(n)) = blocks_[n] * v.segment<3>(3 * Index(n));
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

// Mirrors legacy APolarSite::HistdU(): per-site relative change of the
// induced dipole between successive iterates, with a small-magnitude
// fallback to avoid dividing by a near-zero dipole (e.g. on the first
// iteration from x = 0). Legacy's threshold is 1e-10 e*nm; it is
// converted to bohr since dipoles here are in atomic units. The
// convergence tolerance epstol (SolveWithJOR) is a ratio and needs no
// conversion.
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
// in the structure of legacy PolarBackground's SOR (see
// EwaldSitePolarizabilityBlocks). For A = D + R:
//   x_new = x_old + omega * D^-1 * (b - A*x_old)
// Unlike CG this does not require A to be positive definite; it
// converges if the spectral radius of (I - omega*D^-1*A) is below 1,
// and it does converge on the small indefinite test case where CG fails.
//
// Convergence criterion as in legacy (step III.D): the per-site relative
// dipole change (SiteRelativeDipoleChange) must be <= epstol = 1e-3 for
// every site, or the average must be < 0.1*epstol. The global relative
// residual ||b-Ax||/||b|| is logged for information only.
//
// match_legacy_first_step = true uses omega = 1 for the first iteration
// and the configured omega afterwards, replicating legacy's unrelaxed
// InduceDirect() followed by wSOR-relaxed Induce() calls. The iterates
// are then directly comparable to legacy's (up to the bohr/nm dipole
// conversion).
JorResult SolveWithJOR(const EwaldPeriodicDipoleOperator& op,
                       const EwaldSitePolarizabilityBlocks& site_p,
                       const Eigen::VectorXd& b, Index max_iter, double omega,
                       Logger& log,
                       std::chrono::steady_clock::time_point t_start,
                       bool match_legacy_first_step) {
  auto elapsed_s = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         t_start)
        .count();
  };

  const double kEpsTol = 1e-3;  // legacy's hardcoded value

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
        << ", informational only) (" << elapsed_s() << "s)" << std::flush;

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
