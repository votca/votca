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

// ===========================================================================
// STATUS: new, in-progress home for Hirshfeld real-space partitioning, the
// weight-function scheme chosen for CDFT charge (and, later, spin)
// constraints -- see the design discussion this branch grew out of for why
// Hirshfeld over Becke/SSW (Becke-family partitions are known to give
// unphysical atomic charges when atoms differ significantly in size,
// confirmed via CP2K's own move away from plain Becke for this purpose;
// Hirshfeld's promolecular-density weighting avoids that).
//
// A Hirshfeld weight for atom i is w_i(r) = rho_i^free(r) / sum_j
// rho_j^free(r), where rho_i^free is atom i's own ISOLATED, spherically-
// averaged reference density (DFTEngine::ComputeHirshfeldReferenceDensities,
// via RunAtomicDFT_unrestricted with use_hunds_rule_occupation=true) -- this
// class never touches that atomic-SCF machinery itself (which stays private
// to DFTEngine); it only ever receives the resulting densities as an
// explicit argument, matching DFTGradient's own established pattern in this
// branch: static methods, explicit arguments only, no reaching into
// DFTEngine's private state.
//
// Deliberately NOT yet implemented: the sum-over-atoms denominator, the
// full AO-basis weight matrix W_i (the projection actually needed for the
// Lagrange-multiplier potential added to the Fock matrix), or anything
// about the constraint/Lagrange-multiplier optimization itself. This
// commit is only the first, single-atom building block: evaluating one
// isolated reference density, at one real-space point, given that atom's
// own (small, atom-only) basis re-centered on its real position in the
// molecule -- no embedding into the full molecule's AO basis at all (that
// question, and why it does not arise here, was worked out separately:
// unlike DFTEngine::AtomicGuess's own SAD-guess construction, Hirshfeld
// never needs a molecule-sized combined object).
// ===========================================================================

#pragma once
#ifndef VOTCA_XTP_HIRSHFELDPARTITION_H
#define VOTCA_XTP_HIRSHFELDPARTITION_H

// Local VOTCA includes
#include "aobasis.h"

namespace votca {
namespace xtp {

class HirshfeldPartition {
 public:
  /// Evaluates one isolated-atom reference density (in its own small,
  /// atom-only basis, already re-centered -- via atom_basis's own Fill()
  /// call -- on that atom's real position in the molecule, NOT wherever
  /// the isolated calculation that produced reference_density happened to
  /// place it) at a single real-space point: rho_i^free(point) =
  /// sum_munu P_munu phi_mu(point) phi_nu(point), the standard AO-basis
  /// density-at-a-point formula (matching the same quadratic-form pattern
  /// already used throughout this branch's own XC grid-integration code,
  /// e.g. PulayGradient/GridWeightGradient's own rho = ao.values.dot(DMAT
  /// * ao.values)).
  ///
  /// Deliberately evaluates every shell in atom_basis directly (via each
  /// AOShell's own EvalAOspace), rather than going through GridBox's
  /// "significant shells" filtering (GridBox::CalcAOValues's own
  /// approach) -- atom_basis is already small (one atom's own basis
  /// functions only), so that optimization, built for a full molecule's
  /// much larger basis, does not help here and would only add complexity.
  static double EvaluateAtomicDensity(const AOBasis& atom_basis,
                                       const Eigen::MatrixXd& reference_density,
                                       const Eigen::Vector3d& point);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_HIRSHFELDPARTITION_H
