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

// Hirshfeld real-space partitioning, used as the weight function for CDFT
// charge constraints. Hirshfeld is preferred over Becke/SSW partitioning,
// which gives unphysical atomic charges for atoms of very different size
// (cf. CP2K).
//
// The weight of atom i is w_i(r) = rho_i^free(r) / sum_j rho_j^free(r),
// where rho_i^free is the spherically averaged density of the isolated
// atom (DFTEngine::ComputeHirshfeldReferenceDensities). The reference
// densities are passed in explicitly; as in DFTGradient, all methods are
// static and do not access DFTEngine internals. Each reference density is
// evaluated in its own single-atom basis centered on the atom, so no
// embedding into the molecular AO basis is needed.

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
  /// Reference of one atom: its single-atom basis centered on the atom's
  /// position in the molecule and the reference density matrix of its
  /// element (identical for all atoms of an element).
  struct AtomicReference {
    AOBasis basis;
    Eigen::MatrixXd density;
  };

  /// Builds one AtomicReference per atom of mol (each needs its own
  /// basis center). Call once per geometry and reuse for all grid points.
  static std::vector<AtomicReference> BuildAtomicReferences(
      const QMMolecule& mol, const std::string& basisset_name,
      const std::map<std::string, Eigen::MatrixXd>& reference_densities);

  /// Isolated-atom reference density at point,
  /// rho_i^free(point) = sum_munu P_munu phi_mu(point) phi_nu(point),
  /// with atom_basis centered on the atom in the molecule.
  ///
  /// All shells are evaluated (AOShell::EvalAOspace) without the
  /// significant-shell screening of GridBox, which does not pay off for
  /// a single-atom basis.
  static double EvaluateAtomicDensity(const AOBasis& atom_basis,
                                      const Eigen::MatrixXd& reference_density,
                                      const Eigen::Vector3d& point);

  /// Spatial gradient grad_r rho_i^free(point) of an isolated-atom
  /// reference density, needed for the weight derivative in the CDFT
  /// forces (see EvaluateWeightGradient):
  /// grad_rho(r) = 2 * derivatives^T * (P * values), using the
  /// derivatives returned by AOShell::EvalAOspace.
  static Eigen::Vector3d EvaluateAtomicDensityGradient(
      const AOBasis& atom_basis, const Eigen::MatrixXd& reference_density,
      const Eigen::Vector3d& point);

  /// Hirshfeld weight w_target(point) = rho_target(point) /
  /// sum_j rho_j(point), with the sum over all atoms. Returns 0 where
  /// the denominator is negligible (far from all atoms).
  static double EvaluateWeight(const std::vector<AtomicReference>& atoms,
                               Index target_atom_index,
                               const Eigen::Vector3d& point);

  /// d w_target(point) / d R_A with A = differentiate_atom_index, for the
  /// CDFT forces. Since rho_j(r - R_j) depends only on R_j,
  /// d rho_j/d R_A = -delta_jA grad_r rho_A, and the quotient rule gives
  ///   d w_target/d R_A = [w_target(point) - 1_{A==target}]
  ///                      * grad_rho_A(point) / rho_tot(point).
  /// This is nonzero for every atom A, since all atoms enter rho_tot.
  /// For a multi-atom fragment, sum over the fragment atoms as
  /// target_atom_index, as BuildCDFTConstraint does for BuildWeightMatrix.
  static Eigen::Vector3d EvaluateWeightGradient(
      const std::vector<AtomicReference>& atoms, Index target_atom_index,
      Index differentiate_atom_index, const Eigen::Vector3d& point);

  /// Gradient grad_r w_target(point) with respect to the evaluation
  /// point (EvaluateWeightGradient differentiates with respect to an
  /// atom position instead). Needed for the grid-point-translation force
  /// term: grid points move rigidly with their owner atom
  /// (r_p = R_owner + const, see PulayGradient), so quantities evaluated
  /// at r_p depend on R_owner through r_p.
  ///
  ///   grad_r w(r) = [grad_r rho_target(r) - w(r) * grad_r rho_tot(r)]
  ///                 / rho_tot(r)
  static Eigen::Vector3d EvaluatePointWeightGradient(
      const std::vector<AtomicReference>& atoms, Index target_atom_index,
      const Eigen::Vector3d& point);

  /// CDFT force term from the derivative of the SSW grid weights,
  /// d(Tr[D*W_c])/dR for all atoms of mol (Natoms x 3). The trace is
  /// contracted inside the grid loop, without forming d(W_c)/dR, as in
  /// Vxc_Potential::GridWeightGradient. Implemented separately from that
  /// function (duplicating the SSW functions) to leave the validated XC
  /// gradient code unchanged.
  ///
  /// density_matrix is the converged molecular density in the full AO
  /// basis, from which rho_molecule(point) is evaluated; w_c(point) comes
  /// from the Hirshfeld references in atoms.
  static Eigen::MatrixXd GridWeightDerivativeContribution(
      const std::vector<AtomicReference>& atoms, Index target_atom_index,
      const Eigen::MatrixXd& density_matrix, const QMMolecule& mol,
      const Vxc_Grid& grid);

  /// Pulay (basis-function derivative) and grid-point-translation terms
  /// of d(Tr[D*W_c])/dR (Natoms x 3), computed in one grid loop. The
  /// translation term arises because grid points move with their owner
  /// atom (see EvaluatePointWeightGradient).
  ///
  /// Follows the bookkeeping and sign convention of PulayGradient,
  /// -weight * potential * temp(mu) * ao.derivatives.row(mu), with
  /// w_c(point) in place of xc.df_drho; kept separate from it for the
  /// same reason as GridWeightDerivativeContribution.
  static Eigen::MatrixXd PulayAndTranslationContribution(
      const std::vector<AtomicReference>& atoms, Index target_atom_index,
      const Eigen::MatrixXd& density_matrix, const AOBasis& full_dftbasis,
      const Vxc_Grid& grid);

  /// Weight-function derivative term of the CDFT forces (Natoms x 3):
  /// sum_p omega_p * [d w_c(r_p)/d R_A] * rho_molecule(r_p), with
  /// omega_p and rho_molecule(r_p) held fixed (product-rule term of
  /// omega_p*w_c(r_p)*rho_molecule(r_p)). rho_molecule is evaluated from
  /// the converged density matrix in the full AO basis.
  static Eigen::MatrixXd WeightFunctionDerivativeContribution(
      const std::vector<AtomicReference>& atoms, Index target_atom_index,
      const Eigen::MatrixXd& density_matrix, const AOBasis& full_dftbasis,
      const Vxc_Grid& grid);

  /// Full d(Tr[D*W_c])/dR for all atoms of mol: the sum of the four
  /// product-rule terms of d[omega_p * w_c(r_p) * rho_molecule(r_p)]/dR
  /// (SSW grid weight, weight function, Pulay, grid-point translation).
  /// This is a gradient, not a force; the force of constraint c is
  /// F_i^c = -lambda_c * dN_c/dR_i (Wu-Van Voorhis Lagrangian).
  static Eigen::MatrixXd ComputeCDFTForceContribution(
      const std::vector<AtomicReference>& atoms, Index target_atom_index,
      const Eigen::MatrixXd& density_matrix, const QMMolecule& mol,
      const AOBasis& full_dftbasis, const Vxc_Grid& grid);

  /// AO-basis weight matrix for the Lagrange-multiplier potential,
  /// W_i,munu = integral w_i(r) phi_mu(r) phi_nu(r) dr, integrated
  /// numerically box by box as in Vxc_Potential::IntegrateVXC. grid
  /// should be the grid of the molecular SCF.
  ///
  /// AO values are skipped at points where EvaluateWeight returns
  /// exactly 0 (performance only).
  static Eigen::MatrixXd BuildWeightMatrix(
      const std::vector<AtomicReference>& atoms, Index target_atom_index,
      const AOBasis& full_dftbasis, const Vxc_Grid& grid);

  /// One CDFT constraint. The potential added to the alpha (beta) Fock
  /// matrix is lambda * spin_alpha_coefficient (spin_beta_coefficient)
  /// * weight_matrix.
  ///
  /// A charge constraint, Tr[(P_alpha + P_beta) * W], uses +1/+1 (the
  /// only kind implemented). A spin constraint, Tr[(P_alpha - P_beta) * W],
  /// would use +1/-1 without changes to the Fock-matrix assembly.
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
