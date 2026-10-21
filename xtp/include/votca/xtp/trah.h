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
#ifndef VOTCA_XTP_TRAH_H
#define VOTCA_XTP_TRAH_H

// Standard includes
#include <functional>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief Trust-region augmented-Hessian (TRAH) step for orbital rotations,
 * after Helmich-Paris, J. Chem. Phys. 154, 164104 (2021).
 *
 * Shared by the closed-shell (ConvergenceAcc) and unrestricted
 * (UKSConvergenceAcc) SCF. Orbital rotations are parametrised as
 * C' = C exp(kappa) with kappa_ia = x_ia (i occupied, a virtual) and
 * kappa_ai = -x_ia; the vector x is flattened as x(i * n_virt + a), with the
 * channels of an unrestricted calculation stacked (alpha first). The caller
 * supplies the energy gradient g = dE/dx, the action of the energy Hessian
 * d^2E/dx dx on a block of vectors, and a positive diagonal estimate of the
 * Hessian for preconditioning.
 *
 * The step minimises the quadratic model g.x + 1/2 x.H.x within |x| <= trust
 * radius. It is the lowest eigenvector (mu, [v0; v]) of the scaled augmented
 * Hessian [[0, a g^T], [a g, H]], x = v / (a v0), which satisfies
 * (H - mu) x = -g. All a are handled in one Davidson subspace: only the first
 * row and column of the projected matrix depend on a, so searching a costs
 * no Hessian products. a = 1 is used when its step lies inside the trust
 * radius; otherwise a is increased until |x| lies in
 * [lower_fraction, 1] x trust radius.
 */
namespace trah {

struct StepOptions {
  // The micro-iterations stop when |(H - mu) x + g| <= eta |g| with the
  // forcing term eta = min(max_forcing, |g|), so the macro-iterations
  // converge quadratically close to the minimum.
  double max_forcing = 0.1;
  double residual_floor = 1e-10;
  Index max_iterations = 50;
  Index max_subspace = 60;
  double alpha_max = 1e4;
  double lower_fraction = 0.8;
};

struct StepResult {
  Eigen::VectorXd x;
  // g.x + 1/2 x.H.x of the step taken
  double predicted = 0.0;
  double alpha = 1.0;
  double mu = 0.0;
  Index iterations = 0;
  Index products = 0;
  double residual = 0.0;
  bool converged = false;
};

/// Hessian applied to the columns of V (n x k), returning n x k.
using HessianProduct = std::function<Eigen::MatrixXd(const Eigen::MatrixXd&)>;

StepResult Step(const Eigen::VectorXd& g, const Eigen::VectorXd& diag,
                const HessianProduct& hessian, double trust_radius,
                const StepOptions& opt = StepOptions());

/// Fletcher's trust-radius update. Returns whether the step is accepted
/// (actual and predicted energy changes of the same sign, or both below
/// noise_floor in magnitude) and adjusts the radius: x0.7 for a rejected
/// step or r <= 0.25, x1.2 (up to max_radius) for r > 0.75, with
/// r = actual / predicted.
class TrustRadius {
 public:
  explicit TrustRadius(double radius = 0.5, double max_radius = 1.0)
      : radius_(radius), max_radius_(max_radius) {}
  double Radius() const { return radius_; }
  double LastRatio() const { return ratio_; }
  bool Update(double predicted, double actual);

  static constexpr double noise_floor = 1e-10;

 private:
  double radius_;
  double max_radius_;
  double ratio_ = 0.0;
};

/// C exp(kappa) for the occupied-virtual rotation x (n_occ x n_virt values,
/// layout above), exact via the singular value decomposition of the
/// occupied-virtual block X = U s V^T:
///   exp(kappa) = [[1 + U (cos s - 1) U^T,  U sin s V^T],
///                 [-V sin s U^T,           1 + V (cos s - 1) V^T]].
/// S-orthonormal orbitals stay S-orthonormal.
Eigen::MatrixXd RotateOrbitals(const Eigen::MatrixXd& C,
                               const Eigen::VectorXd& x, Index nocc);

/// Occupied-virtual block of a matrix in the MO basis, flattened as
/// out(i * n_virt + a) = M(i, nocc + a).
Eigen::VectorXd OccVirt(const Eigen::MatrixXd& M, Index nocc);

}  // namespace trah
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TRAH_H
