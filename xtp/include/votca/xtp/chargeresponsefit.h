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
#ifndef VOTCA_XTP_CHARGERESPONSEFIT_H
#define VOTCA_XTP_CHARGERESPONSEFIT_H

// Standard includes
#include <array>
#include <string>
#include <vector>

// Local VOTCA includes
#include "aobasis.h"
#include "eigen.h"
#include "hirshfeldpartition.h"
#include "vxc_grid.h"

namespace votca {
namespace tools {
class Property;
}
namespace xtp {
class Logger;
class Orbitals;

/**
 * \brief Parametrization of a charge-response molecule (see ChargeResponse)
 * from the static density response of a DFT ground state.
 *
 * The perturbations are the Hirshfeld-localized potentials w_k (one per
 * atom; sum_k w_k = 1) and three uniform fields. With the AO matrices W_k of
 * the weights and the first-order density matrices dP of the perturbations
 * (electron potential energy +w_k, resp. +r_c, per unit strength):
 *
 *   Pi_jk    = -Tr(W_j dP_k)                     (dq = -Pi phi)
 *   alpha_ac = -Tr(r_a dP_c)                     (molecular polarizability)
 *   alpha^res = alpha - X^T Pi X                 (X: positions - centre)
 *
 * The residual polarizability is distributed over the atoms as one tensor
 * scaled by the traces of the Hirshfeld atomic polarizabilities
 * alpha_k,ac = -int w_k (r - R_k)_a drho_c. Fit obtains dP by central
 * finite differences of xtpdft SCFs; the partition integrals and the
 * analysis are separate for testing.
 */
class ChargeResponseFit {
 public:
  struct Result {
    Eigen::MatrixXd pi;         // symmetrized, charge conserving
    double pi_asymmetry = 0.0;  // ||Pi - Pi^T|| / ||Pi|| before
    Eigen::VectorXd pi_eigenvalues;
    Eigen::Matrix3d alpha;                         // molecular, a.u.
    Eigen::Matrix3d alpha_ct;                      // X^T Pi X
    Eigen::Matrix3d alpha_res;                     // alpha - alpha_ct
    Eigen::VectorXd res_weights;                   // sum 1
    std::vector<Eigen::Matrix3d> alpha_res_atoms;  // weights * alpha_res
    std::vector<Eigen::Matrix3d> alpha_atoms;      // Hirshfeld, a.u.
  };

  struct FitOptions {
    /// coupled-perturbed KS (DFTEngine::StaticResponse); otherwise central
    /// finite differences of SCFs
    bool analytic = true;
    /// relative CPKS residual
    double tolerance = 1e-7;
    /// finite differences: strength of each perturbation
    double step = 1e-3;
    /// with an inversion centre (FindInversion) only one atom of each pair
    /// is perturbed
    bool use_inversion = true;
  };

  /// The fit with xtpdft (dftoptions: a dftpackage tree with name xtp):
  /// the ground state of orb (atoms, charge and spin set; closed shell),
  /// which orb holds on return, then the first-order densities of all
  /// perturbations, analytically or as dP = (P(+step) - P(-step)) / (2
  /// step) from SCFs warm-started from the ground state. The Hirshfeld
  /// references are xtpdft's (DFTEngine::HirshfeldReferences).
  static Result Fit(const tools::Property& dftoptions, Orbitals& orb,
                    const FitOptions& fit, Logger& log);

  /// Inversion through the centroid of the atoms: atom_map[k] is the image
  /// of atom k, and T the AO matrix of the operation (a signed permutation,
  /// (-1)^l per shell), so that the response to w_k gives the one to
  /// w_atom_map[k] as T dP T^T. atom_map is empty if an atom has no image
  /// of the same element within tol (bohr).
  struct Inversion {
    std::vector<Index> atom_map;
    Eigen::MatrixXd T;
  };
  static Inversion FindInversion(const QMMolecule& mol, const AOBasis& basis,
                                 double tol = 1e-4);

  /// The Hirshfeld weight matrices W_k = int w_k chi chi of all atoms, in
  /// one pass over the grid (the same integrals as
  /// HirshfeldPartition::BuildWeightMatrix).
  static std::vector<Eigen::MatrixXd> WeightMatrices(
      const std::vector<HirshfeldPartition::AtomicReference>& atoms,
      const AOBasis& basis, const Vxc_Grid& grid);

  /// Hirshfeld atomic polarizabilities alpha_k,ac = -int w_k (r - R_k)_a
  /// drho_c(r), with drho_c the density of dP_field[c] (per unit field).
  static std::vector<Eigen::Matrix3d> AtomicPolarizabilities(
      const std::vector<HirshfeldPartition::AtomicReference>& atoms,
      const std::vector<Eigen::Vector3d>& positions, const AOBasis& basis,
      const Vxc_Grid& grid, const std::array<Eigen::MatrixXd, 3>& dP_field);

  /// Pi, alpha and the residual from the responses: dP_potential[k] to
  /// +w_k, dP_field[c] to +r_c, dipole[a] the AO matrices of r_a (any
  /// origin; the responses conserve the electron number). Throws if Pi is
  /// not positive semidefinite or alpha^res not positive definite.
  static Result Analyse(const std::vector<Eigen::Vector3d>& positions,
                        const std::vector<Eigen::MatrixXd>& W,
                        const std::vector<Eigen::MatrixXd>& dP_potential,
                        const std::array<Eigen::MatrixXd, 3>& dipole,
                        const std::array<Eigen::MatrixXd, 3>& dP_field,
                        const std::vector<Eigen::Matrix3d>& alpha_atoms);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CHARGERESPONSEFIT_H
