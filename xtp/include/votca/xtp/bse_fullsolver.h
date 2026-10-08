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
#ifndef VOTCA_XTP_BSE_FULLSOLVER_H
#define VOTCA_XTP_BSE_FULLSOLVER_H

// Standard includes
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

// Third party includes
#include <boost/format.hpp>

// Local VOTCA includes
#include "eigen.h"
#include "logger.h"
#include <votca/tools/eigensystem.h>

namespace votca {
namespace xtp {

/// Scales each column pair of (X, Y) to X^T X - Y^T Y = 1. Columns whose
/// norm is not positive (the zero vectors of unconverged roots, or a root of
/// the wrong sign) are not scaled to nan: zero ones stay zero, negative ones
/// are scaled with |X^T X - Y^T Y|. Returns the number of such columns.
inline Index NormalizeExcitationVectors(Eigen::MatrixXd& X,
                                        Eigen::MatrixXd& Y) {
  Index problems = 0;
  for (Index i = 0; i < X.cols(); ++i) {
    const double norm = X.col(i).squaredNorm() - Y.col(i).squaredNorm();
    const double scale = X.col(i).squaredNorm() + Y.col(i).squaredNorm();
    if (norm > 1e-12 * scale && scale > 0.0) {
      const double f = 1.0 / std::sqrt(norm);
      X.col(i) *= f;
      Y.col(i) *= f;
      continue;
    }
    ++problems;
    if (std::abs(norm) > 1e-12 * scale && scale > 0.0) {
      const double f = 1.0 / std::sqrt(std::abs(norm));
      X.col(i) *= f;
      Y.col(i) *= f;
    }
  }
  return problems;
}

/**
 * \brief Lowest roots of the full (non-TDA) BSE
 *
 *   [ A  B ] [X]       [ 1  0 ] [X]
 *   [ B  A ] [Y] = w   [ 0 -1 ] [Y]
 *
 * for real symmetric A, B with A-B and A+B positive definite, in the
 * equivalent Hermitian form of Stratmann, Scuseria and Frisch
 * (J. Chem. Phys. 109, 8218 (1998)):
 *
 *   (A+B)|X+Y> = w |X-Y>,   (A-B)|X-Y> = w |X+Y>.
 *
 * One subspace V carries both |X+Y> and |X-Y>. Per iteration, A and B are
 * applied once each to the new vectors; the projected problem
 * M- M+ x = w^2 x (M+- = V^T (A+-B) V) is solved through the Cholesky factor
 * of M-, which is symmetric and gives real roots. Restarts collapse the
 * subspace onto the current Ritz vectors without new operator products.
 *
 * Throws std::runtime_error if the projected A-B is not positive definite or
 * a projected root is not positive (an instability of the reference); callers
 * can then fall back to a general non-Hermitian solver.
 */
class FullBSEDavidson {
 public:
  struct Options {
    double tolerance = 1e-4;  // residual norm (Ha) per normalised vector
    Index max_iterations = 50;
    Index max_subspace = 0;  // 0: 10 * nroots
  };

  FullBSEDavidson(Logger& log, Options opt) : log_(log), opt_(opt) {}

  Index Iterations() const { return iterations_; }
  Index OperatorApplications() const { return applications_; }

  template <class OpA, class OpB>
  tools::EigenSystem Solve(const OpA& A, const OpB& B, Index nroots) {
    const auto start = std::chrono::steady_clock::now();
    const Index n = A.rows();
    nroots = std::min(nroots, n);
    // A few guard roots above the requested ones, converged to a looser
    // tolerance: without them the highest requested root can converge to a
    // state above a root the subspace has not seen yet.
    const Index nwork = std::min(n, nroots + std::max<Index>(2, nroots / 4));
    const Index max_subspace = std::min(
        n, std::max(opt_.max_subspace > 0 ? opt_.max_subspace : 10 * nroots,
                    6 * nwork));

    // Diagonal estimate w_i^2 = (a-b)(a+b) for the initial vectors and the
    // preconditioner.
    const Eigen::VectorXd adiag = A.diagonal();
    const Eigen::VectorXd bdiag = B.diagonal();
    std::vector<Index> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    auto estimate = [&](Index i) {
      return std::max(0.0, (adiag(i) - bdiag(i)) * (adiag(i) + bdiag(i)));
    };
    std::stable_sort(order.begin(), order.end(), [&](Index i, Index j) {
      return estimate(i) < estimate(j);
    });
    const Index nguess = std::min(n, 2 * nwork);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, nguess);
    for (Index i = 0; i < nguess; ++i) {
      V(order[std::size_t(i)], i) = 1.0;
    }
    Eigen::MatrixXd P;  // (A+B) V
    Eigen::MatrixXd Q;  // (A-B) V
    Extend(A, B, V, P, Q, V);

    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Full BSE (Hermitian form, " << n << " pairs, "
        << nroots << " roots, tolerance " << opt_.tolerance << ")"
        << std::flush;

    Eigen::VectorXd omega;
    Eigen::MatrixXd x;  // subspace coefficients of |X+Y>
    Eigen::MatrixXd y;  // subspace coefficients of |X-Y>
    Eigen::MatrixXd XpY;
    Eigen::MatrixXd XmY;
    bool converged = false;
    for (iterations_ = 1; iterations_ <= opt_.max_iterations; ++iterations_) {
      ProjectedRoots(V, P, Q, nwork, omega, x, y);
      XpY = V * x;
      XmY = V * y;
      const Eigen::MatrixXd R1 = P * x - XmY * omega.asDiagonal();
      const Eigen::MatrixXd R2 = Q * y - XpY * omega.asDiagonal();
      const Index nr = omega.size();
      Eigen::VectorXd res(nr);
      for (Index i = 0; i < nr; ++i) {
        res(i) = std::max(R1.col(i).norm() / XpY.col(i).norm(),
                          R2.col(i).norm() / XmY.col(i).norm());
      }
      // Requested roots to the tolerance. Guard roots are refined towards
      // 10 x tolerance, but once the requested roots are converged a guard
      // within 100 x tolerance is enough to stop: a guard in a degenerate
      // set cut by the guard window can otherwise take several more
      // iterations without changing the requested roots.
      auto Done = [&](Index i) {
        return res(i) < (i < nroots ? 1.0 : 10.0) * opt_.tolerance;
      };
      Index nconv = 0;
      for (Index i = 0; i < nroots; ++i) {
        nconv += Done(i) ? 1 : 0;
      }
      const double guard_res =
          (nr > nroots) ? res.tail(nr - nroots).maxCoeff() : 0.0;
      const bool all_done =
          nconv == nroots && guard_res < 100.0 * opt_.tolerance;
      XTP_LOG(Log::error, log_)
          << TimeStamp()
          << boost::format(
                 " iter %1$3d  subspace %2$5d  max residual "
                 "%3$4.2e  converged %4$4d/%5$d  guard roots %6$4.2e") %
                 iterations_ % V.cols() % res.head(nroots).maxCoeff() % nconv %
                 nroots % guard_res
          << std::flush;
      if (all_done) {
        converged = true;
        break;
      }

      // Corrections from the diagonal approximation A ~ a, B ~ b of
      //   (A+B) dp - w dm = -r1,  (A-B) dm - w dp = -r2,
      // solved per pair (2x2 system with determinant a^2 - b^2 - w^2).
      std::vector<Eigen::VectorXd> directions;
      for (Index i = 0; i < nr; ++i) {
        if (Done(i)) {
          continue;
        }
        const double w = omega(i);
        Eigen::ArrayXd det =
            adiag.array().square() - bdiag.array().square() - w * w;
        det = det.sign() * det.abs().max(1e-4) +
              (det == 0.0).cast<double>() * 1e-4;
        const Eigen::ArrayXd r1 = R1.col(i).array();
        const Eigen::ArrayXd r2 = R2.col(i).array();
        directions.push_back(
            (-((adiag - bdiag).array() * r1 + w * r2) / det).matrix());
        directions.push_back(
            (-(w * r1 + (adiag + bdiag).array() * r2) / det).matrix());
      }
      Eigen::MatrixXd T(n, Index(directions.size()));
      for (Index i = 0; i < T.cols(); ++i) {
        T.col(i) = directions[std::size_t(i)];
      }

      if (V.cols() + T.cols() > max_subspace) {
        // thick restart: keep the Ritz pairs of twice as many roots
        Eigen::VectorXd omega_keep;
        Eigen::MatrixXd x_keep;
        Eigen::MatrixXd y_keep;
        ProjectedRoots(V, P, Q, std::min(2 * nwork, V.cols()), omega_keep,
                       x_keep, y_keep);
        Collapse(V, P, Q, x_keep, y_keep);
      }
      const Index before = V.cols();
      Extend(A, B, T, P, Q, V);
      if (V.cols() == before) {
        XTP_LOG(Log::error, log_)
            << TimeStamp()
            << " Full BSE: no new search directions, stopping with max "
               "residual "
            << res.head(nroots).maxCoeff() << std::flush;
        break;
      }
    }
    if (!converged) {
      ProjectedRoots(V, P, Q, nwork, omega, x, y);
      XpY = V * x;
      XmY = V * y;
      XTP_LOG(Log::error, log_)
          << TimeStamp() << " WARNING: full BSE not converged after "
          << opt_.max_iterations << " iterations" << std::flush;
    }

    // X = (|X+Y> + |X-Y>)/2, Y = (|X+Y> - |X-Y>)/2; x^T y = 1 gives
    // X^T X - Y^T Y = 1.
    tools::EigenSystem result;
    const Index nout = std::min(nroots, Index(omega.size()));
    result.eigenvalues() = omega.head(nout);
    result.eigenvectors() = 0.5 * (XpY + XmY).leftCols(nout);
    result.eigenvectors2() = 0.5 * (XpY - XmY).leftCols(nout);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " Full BSE done: " << std::min(iterations_, opt_.max_iterations)
        << " iterations, " << applications_ << " vectors through A and B, "
        << seconds << " s" << std::flush;
    return result;
  }

 private:
  // Orthonormalises T against V and itself (two Gram-Schmidt passes), drops
  // directions without a significant new component, applies A and B to the
  // rest and appends them to V, P and Q.
  template <class OpA, class OpB>
  void Extend(const OpA& A, const OpB& B, Eigen::MatrixXd T, Eigen::MatrixXd& P,
              Eigen::MatrixXd& Q, Eigen::MatrixXd& V) {
    const bool fresh = (P.cols() == 0);
    Eigen::MatrixXd basis = fresh ? Eigen::MatrixXd(V.rows(), 0) : V;
    std::vector<Index> kept;
    Eigen::MatrixXd accepted(V.rows(), 0);
    for (Index i = 0; i < T.cols(); ++i) {
      Eigen::VectorXd t = T.col(i);
      const double norm0 = t.norm();
      if (norm0 == 0.0) {
        continue;
      }
      for (int pass = 0; pass < 2; ++pass) {
        if (basis.cols() > 0) {
          t -= basis * (basis.transpose() * t);
        }
        if (accepted.cols() > 0) {
          t -= accepted * (accepted.transpose() * t);
        }
      }
      const double norm = t.norm();
      if (norm < 1e-8 * norm0 || norm < 1e-12) {
        continue;
      }
      accepted.conservativeResize(Eigen::NoChange, accepted.cols() + 1);
      accepted.col(accepted.cols() - 1) = t / norm;
    }
    if (accepted.cols() == 0) {
      if (fresh) {
        V.resize(V.rows(), 0);
      }
      return;
    }
    const Eigen::MatrixXd At = A.matmul(accepted);
    const Eigen::MatrixXd Bt = B.matmul(accepted);
    applications_ += accepted.cols();
    if (fresh) {
      V = accepted;
      P = At + Bt;
      Q = At - Bt;
      return;
    }
    const Index m = V.cols();
    const Index k = accepted.cols();
    V.conservativeResize(Eigen::NoChange, m + k);
    P.conservativeResize(Eigen::NoChange, m + k);
    Q.conservativeResize(Eigen::NoChange, m + k);
    V.rightCols(k) = accepted;
    P.rightCols(k) = At + Bt;
    Q.rightCols(k) = At - Bt;
  }

  void ProjectedRoots(const Eigen::MatrixXd& V, const Eigen::MatrixXd& P,
                      const Eigen::MatrixXd& Q, Index nroots,
                      Eigen::VectorXd& omega, Eigen::MatrixXd& x,
                      Eigen::MatrixXd& y) const {
    Eigen::MatrixXd Mp = V.transpose() * P;
    Eigen::MatrixXd Mm = V.transpose() * Q;
    Mp = 0.5 * (Mp + Mp.transpose()).eval();
    Mm = 0.5 * (Mm + Mm.transpose()).eval();
    Eigen::LLT<Eigen::MatrixXd> llt(Mm);
    if (llt.info() != Eigen::Success) {
      throw std::runtime_error(
          "FullBSEDavidson: projected A-B is not positive definite");
    }
    const Eigen::MatrixXd L = llt.matrixL();
    Eigen::MatrixXd H = L.transpose() * Mp * L;
    H = 0.5 * (H + H.transpose()).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    const Index nr = std::min(nroots, Index(es.eigenvalues().size()));
    if (es.eigenvalues()(0) <= 0.0) {
      throw std::runtime_error(
          "FullBSEDavidson: non-positive root, A+B is not positive definite "
          "(instability of the reference)");
    }
    omega = es.eigenvalues().head(nr).cwiseSqrt();
    const Eigen::MatrixXd u = es.eigenvectors().leftCols(nr);
    // x = L u / sqrt(w), y = L^-T u sqrt(w): M+ x = w y, M- y = w x, x^T y = 1
    x = L * u * omega.cwiseSqrt().cwiseInverse().asDiagonal();
    y = L.transpose().triangularView<Eigen::Upper>().solve(u) *
        omega.cwiseSqrt().asDiagonal();
  }

  // Replaces the subspace by an orthonormal basis of the current |X+Y> and
  // |X-Y> Ritz vectors; P and Q follow by the same linear combination.
  static void Collapse(Eigen::MatrixXd& V, Eigen::MatrixXd& P,
                       Eigen::MatrixXd& Q, const Eigen::MatrixXd& x,
                       const Eigen::MatrixXd& y) {
    Eigen::MatrixXd S(x.rows(), x.cols() + y.cols());
    S << x, y;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(S);
    const Index rank = qr.rank();
    const Eigen::MatrixXd C = Eigen::MatrixXd(qr.householderQ()).leftCols(rank);
    V = (V * C).eval();
    P = (P * C).eval();
    Q = (Q * C).eval();
  }

  Logger& log_;
  Options opt_;
  Index iterations_ = 0;
  Index applications_ = 0;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BSE_FULLSOLVER_H
