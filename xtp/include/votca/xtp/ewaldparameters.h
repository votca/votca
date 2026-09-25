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
#ifndef VOTCA_XTP_EWALDPARAMETERS_H
#define VOTCA_XTP_EWALDPARAMETERS_H

// Standard includes
#include <string>

// Local VOTCA includes
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/eigen.h"
#include "votca/xtp/ewaldshapecorrection.h"

namespace votca {
namespace xtp {

/**
 * \brief The Ewald convergence parameters a background was converged with.
 *
 * Written into the background's own checkpoint alongside the converged
 * dipoles, and read back by whatever later embeds a region in that
 * background (MM/MM, QM/MM).
 *
 * The point is that these are NOT re-specified by the consuming job.
 * A foreground embedded in a periodic background is only meaningful if
 * it uses the same alpha, the same k-space cutoff and the same
 * summation shape the background was converged with: alpha in
 * particular partitions the interaction between the real- and
 * reciprocal-space sums, so two different values do not describe the
 * same physics at all, and the mismatch would not announce itself --
 * both runs would converge happily to different answers. Carrying them
 * in the checkpoint makes the inconsistency unrepresentable rather than
 * merely discouraged.
 *
 * The box is included for the same reason: the background's lattice is
 * a property of the converged state, not something a later job should
 * re-derive and hope matches.
 */
struct EwaldParameters {
  double alpha = 0.0;  // bohr^-1
  double k_max = 0.0;  // bohr^-1
  double r_min = 0.0;  // bohr
  double field_tol = 0.0;
  double thole_a = 0.0;
  double screening_factor = 0.0;  // dimensionless
  EwaldShape shape = EwaldShape::Cube;
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();

  void WriteToCpt(CheckpointWriter& w) const;
  void ReadFromCpt(CheckpointReader& r);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDPARAMETERS_H
