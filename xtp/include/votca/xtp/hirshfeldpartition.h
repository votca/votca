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

// Standard includes
#include <map>
#include <string>
#include <vector>

// Local VOTCA includes
#include "aobasis.h"
#include "qmmolecule.h"
#include "vxc_grid.h"

namespace votca {
namespace xtp {

class HirshfeldPartition {
 public:
  /// One real atom's own contribution to the Hirshfeld partition: its
  /// small, atom-only basis, re-centered (via its own Fill() call) on
  /// that atom's REAL position in the molecule, together with its
  /// element's reference density matrix (the same matrix is shared,
  /// by value, across every atom of the same element -- only the
  /// basis's own center differs per atom).
  struct AtomicReference {
    AOBasis basis;
    Eigen::MatrixXd density;
  };

  /// Builds one AtomicReference per real atom in mol (NOT one per
  /// unique element -- every atom needs its own, distinctly-centered
  /// basis, even though atoms of the same element share the same
  /// reference density matrix). Intended to be called ONCE per
  /// molecule/SCF, then reused across every grid point EvaluateWeight
  /// is called for -- rebuilding an AOBasis per point would be
  /// wasteful, since only the evaluation point changes from one call
  /// to the next, not the geometry.
  static std::vector<AtomicReference> BuildAtomicReferences(
      const QMMolecule& mol, const std::string& basisset_name,
      const std::map<std::string, Eigen::MatrixXd>& reference_densities);

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

  /// The actual Hirshfeld weight: w_target(point) = rho_target(point) /
  /// sum_j rho_j(point), summing rho_j over EVERY atom in atoms (not
  /// just target_atom_index -- the denominator needs every atom's
  /// contribution, per the Hirshfeld definition). atoms is intended to
  /// be built once via BuildAtomicReferences and reused across many
  /// calls at different points. Returns 0 (rather than dividing by a
  /// near-zero denominator) at points far from every atom, where every
  /// rho_j(point) is negligible -- matching the same
  /// negligible-denominator guard pattern already used for the SSW
  /// grid weights in GridWeightGradient (kNegligibleWOwner there).
  static double EvaluateWeight(const std::vector<AtomicReference>& atoms,
                                Index target_atom_index,
                                const Eigen::Vector3d& point);

  /// The actual AO-basis operator matrix the Lagrange-multiplier
  /// potential needs: W_i,munu = integral w_i(r) phi_mu(r) phi_nu(r) dr,
  /// numerically integrated over full_dftbasis's own molecule-wide
  /// grid -- the same grid, and the same box-by-box
  /// CalcAOValues/AddtoBigMatrix pattern, already used throughout this
  /// branch's own XC potential/Pulay-gradient code
  /// (vxc_potential.cc's IntegrateVXC/PulayGradient). atoms should be
  /// built once via BuildAtomicReferences and passed in here; grid
  /// should be the SAME grid used for the real molecule's own SCF (not
  /// a separately-built one), so this weight matrix is evaluated at
  /// exactly the same points already being used for everything else.
  ///
  /// Skips evaluating AO values entirely at any point where
  /// EvaluateWeight itself returns exactly 0 (points far from every
  /// atom, or on the far side of the molecule from target_atom_index)
  /// -- purely a performance guard, not a correctness one, since a
  /// zero weight contributes nothing to the integral regardless of
  /// what the AO values happen to be there.
  static Eigen::MatrixXd BuildWeightMatrix(
      const std::vector<AtomicReference>& atoms, Index target_atom_index,
      const AOBasis& full_dftbasis, const Vxc_Grid& grid);

  /// One active CDFT constraint: the Lagrange-multiplier potential term
  /// added to the Fock matrix(es) is lambda * spin_alpha_coefficient *
  /// weight_matrix for the alpha channel, and lambda *
  /// spin_beta_coefficient * weight_matrix for the beta channel.
  ///
  /// For a CHARGE constraint (the only kind implemented so far):
  /// spin_alpha_coefficient = spin_beta_coefficient = +1.0, since the
  /// constrained quantity is Tr[(P_alpha + P_beta) * W] -- the SAME
  /// potential added to both spin channels. A SPIN constraint (not yet
  /// implemented, but the reason these two coefficients are stored
  /// separately rather than hardcoding "add to both channels equally")
  /// would instead need +1.0/-1.0, since it constrains
  /// Tr[(P_alpha - P_beta) * W] -- per the design discussion this
  /// struct grew out of, adding spin constraints later should only
  /// ever require constructing a Constraint with different
  /// coefficients, never touching the Fock-matrix-assembly code that
  /// consumes this struct at all.
  struct Constraint {
    Eigen::MatrixXd weight_matrix;
    double target_population;
    double lambda = 0.0;
    double spin_alpha_coefficient = 1.0;
    double spin_beta_coefficient = 1.0;
  };
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_HIRSHFELDPARTITION_H
