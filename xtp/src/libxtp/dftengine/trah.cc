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

// Standard includes
#include <algorithm>
#include <cmath>

// Local VOTCA includes
#include "votca/xtp/trah.h"

namespace votca {
namespace xtp {
namespace trah {

namespace {

// Lowest eigenpair of the projected augmented Hessian
// [[0, a gs^T], [a gs, Hs]]; xs = y / (a y0). False if y0 vanishes (g
// orthogonal to the lowest eigenvector), where no step can be formed.
bool SolveProjected(const Eigen::MatrixXd& Hs, const Eigen::VectorXd& gs,
                    double a, Eigen::VectorXd& xs, double& mu) {
  const Index m = Hs.rows();
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(m + 1, m + 1);
  M.block(0, 1, 1, m) = a * gs.transpose();
  M.block(1, 0, m, 1) = a * gs;
  M.bottomRightCorner(m, m) = Hs;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
  const Eigen::VectorXd y = es.eigenvectors().col(0);
  if (std::abs(y(0)) < 1e-14) {
    return false;
  }
  xs = y.tail(m) / (a * y(0));
  mu = es.eigenvalues()(0);
  return true;
}

// Orthonormalises v against the columns of B (twice) and normalises it.
// False if nothing new is left.
bool Orthonormalize(const Eigen::MatrixXd& B, Eigen::VectorXd& v) {
  const double norm0 = v.norm();
  if (!(norm0 > 0.0)) {
    return false;
  }
  for (int pass = 0; pass < 2; ++pass) {
    if (B.cols() > 0) {
      v -= B * (B.transpose() * v);
    }
  }
  const double norm = v.norm();
  if (norm < 1e-10 * norm0) {
    return false;
  }
  v /= norm;
  return true;
}
}  // namespace

StepResult Step(const Eigen::VectorXd& g, const Eigen::VectorXd& diag,
                const HessianProduct& hessian, double trust_radius,
                const StepOptions& opt, Subspace* subspace) {
  const Index n = g.size();
  StepResult result;
  result.x = Eigen::VectorXd::Zero(n);
  const double gnorm = g.norm();
  if (!(gnorm > 0.0)) {
    result.converged = true;
    return result;
  }
  const double eta = std::min(opt.max_forcing, gnorm);
  const double tolerance = std::max(eta * gnorm, opt.residual_floor);

  Eigen::MatrixXd B(n, 0);
  Eigen::MatrixXd HB(n, 0);
  auto add = [&](Eigen::VectorXd v) {
    if (!Orthonormalize(B, v)) {
      return false;
    }
    const Eigen::MatrixXd Hv = hessian(v);
    ++result.products;
    B.conservativeResize(n, B.cols() + 1);
    B.col(B.cols() - 1) = v;
    HB.conservativeResize(n, HB.cols() + 1);
    HB.col(HB.cols() - 1) = Hv.col(0);
    return true;
  };
  if (subspace != nullptr && subspace->B.cols() > 0 &&
      subspace->B.rows() == n && subspace->g.size() == n &&
      (subspace->g - g).norm() <= 1e-12 * gnorm) {
    B = subspace->B;
    HB = subspace->HB;
    result.reused = true;
  } else {
    add(g / gnorm);
  }

  Eigen::VectorXd xs;
  Eigen::VectorXd gs;
  Eigen::MatrixXd Hs;
  double mu = 0.0;
  double alpha = 1.0;
  for (Index it = 0; it < opt.max_iterations; ++it) {
    result.iterations = it + 1;
    Hs = B.transpose() * HB;
    Hs = 0.5 * (Hs + Hs.transpose()).eval();
    gs = B.transpose() * g;

    // a = 1 if its step is inside the trust radius; otherwise bisection in
    // log(a), with no further Hessian products.
    alpha = 1.0;
    if (!SolveProjected(Hs, gs, alpha, xs, mu)) {
      break;
    }
    if (xs.norm() > trust_radius) {
      double lo = 1.0;
      double hi = opt.alpha_max;
      Eigen::VectorXd xs_hi;
      double mu_hi = 0.0;
      if (SolveProjected(Hs, gs, hi, xs_hi, mu_hi) &&
          xs_hi.norm() > opt.lower_fraction * trust_radius) {
        // even the largest a gives a long step; it is scaled back below
        alpha = hi;
        xs = xs_hi;
        mu = mu_hi;
      } else {
        for (int k = 0; k < 100; ++k) {
          const double a = std::sqrt(lo * hi);
          Eigen::VectorXd xs_a;
          double mu_a = 0.0;
          if (!SolveProjected(Hs, gs, a, xs_a, mu_a)) {
            lo = a;
            continue;
          }
          alpha = a;
          xs = xs_a;
          mu = mu_a;
          const double norm = xs_a.norm();
          if (norm > trust_radius) {
            lo = a;
          } else if (norm < opt.lower_fraction * trust_radius) {
            hi = a;
          } else {
            break;
          }
        }
      }
    }

    // residual of (H - mu) x = -g in the full space
    const Eigen::VectorXd r = HB * xs + g - mu * (B * xs);
    result.residual = r.norm();
    if (result.residual <= tolerance) {
      result.converged = true;
      break;
    }
    if (it + 1 == opt.max_iterations) {
      break;
    }
    if (B.cols() >= opt.max_subspace) {
      // restart on span{x, g}: both are combinations of B, so are their
      // Hessian products
      Eigen::MatrixXd W(B.cols(), 2);
      W.col(0) = xs;
      W.col(1) = gs;
      Eigen::HouseholderQR<Eigen::MatrixXd> qr(W);
      const Eigen::MatrixXd Q =
          qr.householderQ() * Eigen::MatrixXd::Identity(B.cols(), 2);
      B = (B * Q).eval();
      HB = (HB * Q).eval();
    }
    // preconditioned residual; diag - mu > 0 for a positive diag, since
    // mu <= 0 (the projected matrix has a zero on its diagonal)
    Eigen::VectorXd t(n);
    for (Index k = 0; k < n; ++k) {
      t(k) = -r(k) / std::max(diag(k) - mu, 1e-4);
    }
    if (!add(t)) {
      break;
    }
  }

  // predicted change g.x + 1/2 x.H.x, exact in the subspace; a step still
  // outside the trust radius is scaled back onto it
  const double norm = xs.norm();
  const double s = (norm > trust_radius) ? trust_radius / norm : 1.0;
  result.predicted = s * gs.dot(xs) + 0.5 * s * s * xs.dot(Hs * xs);
  result.linear = s * gs.dot(xs);
  result.x = s * (B * xs);
  if (subspace != nullptr) {
    subspace->B = B;
    subspace->HB = HB;
    subspace->g = g;
  }
  result.alpha = alpha;
  result.mu = mu;
  return result;
}

bool TrustRadius::Update(double predicted, double actual, double linear,
                         double step_norm) {
  const double radius_before = radius_;
  const bool accepted = Update(predicted, actual);
  if (!accepted && step_norm > 0.0) {
    const double c = actual - linear;
    double factor = 0.7;
    if (linear < 0.0 && c > 0.0) {
      factor = std::clamp(-linear / (2.0 * c), 0.1, 0.7);
    }
    radius_ = factor * std::min(radius_before, step_norm);
  }
  return accepted;
}

bool TrustRadius::Update(double predicted, double actual) {
  if (std::abs(predicted) < noise_floor && std::abs(actual) < noise_floor) {
    ratio_ = 1.0;
    return true;
  }
  ratio_ = (std::abs(predicted) > 0.0) ? actual / predicted : -1.0;
  if (ratio_ < 0.0) {
    radius_ *= 0.7;
    return false;
  }
  if (ratio_ <= 0.25) {
    radius_ *= 0.7;
  } else if (ratio_ > 0.75) {
    radius_ = std::min(1.2 * radius_, max_radius_);
  }
  return true;
}

Eigen::MatrixXd RotateOrbitals(const Eigen::MatrixXd& C,
                               const Eigen::VectorXd& x, Index nocc) {
  const Index nmo = C.cols();
  const Index nvirt = nmo - nocc;
  if (nocc == 0 || nvirt == 0) {
    return C;
  }
  const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                       Eigen::RowMajor>>
      X(x.data(), nocc, nvirt);
  Eigen::BDCSVD<Eigen::MatrixXd> svd(X,
                                     Eigen::ComputeThinU | Eigen::ComputeThinV);
  const Eigen::MatrixXd& U = svd.matrixU();
  const Eigen::MatrixXd& V = svd.matrixV();
  const Eigen::ArrayXd s = svd.singularValues().array();
  const Eigen::VectorXd cos_m1 = (s.cos() - 1.0).matrix();
  const Eigen::VectorXd sin_s = s.sin().matrix();

  const Eigen::MatrixXd E_oo = Eigen::MatrixXd::Identity(nocc, nocc) +
                               U * cos_m1.asDiagonal() * U.transpose();
  const Eigen::MatrixXd E_ov = U * sin_s.asDiagonal() * V.transpose();
  const Eigen::MatrixXd E_vv = Eigen::MatrixXd::Identity(nvirt, nvirt) +
                               V * cos_m1.asDiagonal() * V.transpose();

  const auto C_occ = C.leftCols(nocc);
  const auto C_virt = C.rightCols(nvirt);
  Eigen::MatrixXd result(C.rows(), nmo);
  result.leftCols(nocc) = C_occ * E_oo - C_virt * E_ov.transpose();
  result.rightCols(nvirt) = C_occ * E_ov + C_virt * E_vv;
  return result;
}

Eigen::VectorXd OccVirt(const Eigen::MatrixXd& M, Index nocc) {
  const Index nvirt = M.cols() - nocc;
  Eigen::VectorXd out(nocc * nvirt);
  for (Index i = 0; i < nocc; ++i) {
    for (Index a = 0; a < nvirt; ++a) {
      out(i * nvirt + a) = M(i, nocc + a);
    }
  }
  return out;
}

}  // namespace trah
}  // namespace xtp
}  // namespace votca
