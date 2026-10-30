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
#ifndef VOTCA_XTP_SCALEDSOLVE_H
#define VOTCA_XTP_SCALEDSOLVE_H

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief Solves A x = b with an Eigen iterative solver for the right-hand
 * side scaled to a largest entry of 1.
 *
 * The solution is the same for a linear system, but the solver sees a
 * right-hand side of order one. Eigen's BiCGSTAB (from version 5) detects a
 * breakdown with thresholds that do not scale with the right-hand side, so
 * small right-hand sides (e.g. times in seconds, about 1e-11, or the
 * residuals of a refinement step) stopped it at once with NumericalIssue
 * and an inaccurate solution. A zero right-hand side gives zero without
 * calling the solver; `solved` tells whether the solver ran (its
 * iterations() and info() refer to this call only then).
 */
template <class IterativeSolver>
Eigen::VectorXd ScaledSolve(const IterativeSolver& solver,
                            const Eigen::VectorXd& b, bool* solved = nullptr) {
  const double scale = b.size() > 0 ? b.cwiseAbs().maxCoeff() : 0.0;
  if (solved != nullptr) {
    *solved = scale > 0.0;
  }
  if (!(scale > 0.0)) {
    return Eigen::VectorXd::Zero(b.size());
  }
  Eigen::VectorXd x = solver.solve(b / scale);
  return x * scale;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SCALEDSOLVE_H
