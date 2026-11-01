/*
 *            Copyright 2009-2024 The VOTCA Development Team
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

// Analytic DFT nuclear gradient terms (nuclear repulsion, RI-J, RI-K),
// assembled from the integral derivatives of libint2_derivative_calls.cc.
// All functions return dE/dR as (Natoms x 3) in Hartree/Bohr.

#pragma once
#ifndef VOTCA_XTP_DFTGRADIENT_H
#define VOTCA_XTP_DFTGRADIENT_H

// Local VOTCA includes
#include "aobasis.h"
#include "qmmolecule.h"

namespace votca {
namespace xtp {

class DFTGradient {
 public:
  /// Derivative of the classical nuclear-nuclear repulsion energy
  /// E_nn = sum over atom pairs A,B with A \< B of Z_A Z_B / R_AB, with
  /// respect to nuclear coordinates.
  /// Returns an (Natoms x 3) matrix in Hartree/Bohr (positions are in
  /// Bohr).
  ///
  /// Uses QMAtom::getNuccharge() (nuclear charge minus ECP core charge,
  /// if an ECP is in use), the effective charge that should enter the
  /// classical repulsion term consistently with an ECP-modified
  /// Hamiltonian.
  static Eigen::MatrixXd NuclearRepulsionDerivative(const QMMolecule& mol);

  /// RI-J (Coulomb) gradient for a fixed density matrix, contracting the
  /// two- and three-centre integral derivatives with the RI fitting
  /// coefficients c = V^-1 d. Returns dE_J/dR (Natoms x 3, Hartree/Bohr).
  ///
  /// The response of c need not be differentiated: E_J(c) = c^T d -
  /// 1/2 c^T V c is stationary in c for any fixed d. The formula is
  /// therefore valid for any density; for forces it must be the converged
  /// SCF density.
  static Eigen::MatrixXd RIJGradient(const Eigen::MatrixXd& density,
                                     const AOBasis& auxbasis,
                                     const AOBasis& dftbasis);

  /// RI-K (exchange) gradient, without the hybrid scaling factor.
  ///
  /// occ_mo_coeffs: (nbf x nocc) occupied MO coefficients C. With
  ///   d_ij(P) = C_i^T (mu nu|P) C_j,  c_ij = V^-1 d_ij,
  ///   E_K = -2 * sum_{i,j} [0.5 * c_ij . d_ij] = -sum_{i,j} c_ij . d_ij
  /// over all ordered pairs of occupied MOs, which equals the exchange
  /// energy of ERIs::CalculateEXX_mos. That uses the symmetric V^-1/2 fit,
  /// which gives the same energy as the V^-1 fit used here
  /// (|V^-1/2 x|^2 = x^T V^-1 x for positive-definite V), so no derivative
  /// of V^-1/2 is needed. As for RIJGradient, c_ij is stationary for any
  /// fixed C, so the formula holds for any C; for forces C must be the
  /// converged occupied MOs.
  static Eigen::MatrixXd RIKGradient(const Eigen::MatrixXd& occ_mo_coeffs,
                                     const AOBasis& auxbasis,
                                     const AOBasis& dftbasis);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTGRADIENT_H
