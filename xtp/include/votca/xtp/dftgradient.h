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

// ===========================================================================
// STATUS: new, in-progress home for DFT gradient (nuclear force) assembly,
// as distinct from the raw integral derivatives in libint2_derivative_calls
// (which this consumes). Currently implements ONLY the nuclear repulsion
// term. Not yet wired into DFTEngine or Orbitals::setForces() -- that is
// a separate, later step once the RI-J (and, if hybrids matter, RI-K) and
// XC gradient terms exist alongside this one.
//
// Planned eventual contents, per the dependency order worked out
// separately: nuclear repulsion (this commit) -> RI-J gradient assembly,
// contracting the already-validated two-/three-center integral
// derivatives with the converged density and fitting coefficients ->
// XC gradient (density-derivative + SSW grid-weight-derivative terms,
// the one piece with no existing shortcut in the codebase) -> total
// gradient, called from DFTEngine once SCF has converged, written to
// Orbitals::setForces(). CDFT's own lambda_i * d(w_i)/dR term rides on
// the XC term's grid-weight-derivative code, once that exists, and is
// deliberately out of scope until the plain-DFT case is complete and
// tested end to end.
// ===========================================================================

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
  /// E_nn = sum_{A<B} Z_A Z_B / R_AB with respect to nuclear coordinates.
  /// Returns an (Natoms x 3) matrix in Hartree/Bohr, consistent with
  /// QMAtom::getPos() being stored in Bohr and with the
  /// Orbitals::setForces() convention -- no unit conversion needed here,
  /// unlike the libint2-derived integral derivatives (which needed a
  /// Bohr/Angstrom fix in their own finite-difference tests): this
  /// function consumes QMAtom positions directly, already in Bohr, and
  /// produces a result already in atomic units end to end.
  ///
  /// Uses QMAtom::getNuccharge() (nuclear charge minus ECP core charge,
  /// if an ECP is in use), the effective charge that should enter the
  /// classical repulsion term consistently with an ECP-modified
  /// Hamiltonian.
  static Eigen::MatrixXd NuclearRepulsionDerivative(const QMMolecule& mol);

  /// RI-J (Coulomb) gradient assembly, contracting the already-validated
  /// two-/three-center integral derivatives with a FIXED density matrix
  /// and the corresponding RI fitting coefficients. Returns (Natoms x 3)
  /// dE_J/dR in Hartree/Bohr, same convention as
  /// NuclearRepulsionDerivative.
  ///
  /// IMPORTANT: "density" here does not need to come from a converged
  /// SCF -- the underlying identity this relies on (that the RI fitting
  /// coefficients c_P are stationary in the fitted energy expression, so
  /// their own response to a geometry change doesn't need to be
  /// differentiated) is a property of the linear least-squares fit
  /// itself (c = V^-1 d makes E_J(c) = c^T d - 1/2 c^T V c stationary in
  /// c by construction, for ANY fixed d), not a consequence of
  /// electronic self-consistency. This function is valid for any fixed
  /// density matrix passed in, converged or not; it is the CALLER's
  /// responsibility to pass in a genuinely converged SCF density when
  /// this is used for a real gradient (as opposed to, e.g., a testing
  /// context where an arbitrary fixed density is a legitimate way to
  /// validate the assembly formula itself).
  static Eigen::MatrixXd RIJGradient(const Eigen::MatrixXd& density,
                                      const AOBasis& auxbasis,
                                      const AOBasis& dftbasis);

  /// RI-K (exchange) gradient assembly.
  ///
  /// HISTORY, worth keeping for the record: a brief intermediate
  /// revision switched this to a "half-transformed" structure (one
  /// index MO, one AO), reasoned via error-prone hand algebra to be
  /// needed to match the real K matrix (see ERIs::CalculateEXX_mos).
  /// That reasoning was WRONG. Settled by directly, numerically
  /// simulating CalculateEXX_mos's actual algorithm (symmetric V^-1/2
  /// fit) against both candidate formulas on several random test
  /// systems: the FULLY-MO-transformed structure below (both indices
  /// occupied MOs) is correct and matches the real energy exactly (to
  /// ~1e-14) once a missing factor of 2 is included -- the
  /// half-transformed version did not match at all, not even up to a
  /// constant factor.
  ///
  /// Separately confirmed (also numerically, to machine precision): for
  /// symmetric positive-definite V, |V^-1/2 x|^2 == x^T V^-1 x EXACTLY,
  /// an algebraic identity -- so the symmetric (production) and
  /// asymmetric (used here) RI fitting schemes give IDENTICAL total
  /// exchange energies, and therefore identical gradients under the
  /// same stationarity argument used throughout this file. No matrix
  /// square root derivative is needed; V^-1 is both simpler AND exact,
  /// once the factor of 2 is included.
  ///
  /// occ_mo_coeffs: (nbf x nocc) occupied MO coefficient matrix. Despite
  /// the name, does NOT need to be genuinely orthonormal SCF occupied
  /// orbitals for this function to be mathematically valid -- same
  /// reasoning as RIJGradient's density argument: the per-(i,j)-pair
  /// fitting coefficients are stationary in the fitted exchange energy
  /// expression for ANY fixed C, a pure linear-algebra property, not a
  /// consequence of SCF self-consistency or orthonormality. (In
  /// PRACTICE, for this to match a real hybrid SCF's exact-exchange
  /// gradient, it must be called with the genuine converged occupied
  /// MOs -- the "doesn't need to be orthonormal" note is about
  /// mathematical validity for testing purposes, per this file's
  /// established pattern, not a license to pass arbitrary MOs in
  /// production use.)
  ///
  /// Energy convention (matches ERIs::CalculateEXX_mos's REAL physical
  /// exchange energy exactly, confirmed numerically, not just up to an
  /// unknown scale):
  ///   d_ij(P) = C_i^T (AO 3-center integral tensor[P]) C_j
  ///   c_ij = V^-1 d_ij   (per (i,j) pair)
  ///   E_K = -2 * sum_{i,j} [0.5 * c_ij . d_ij] = -sum_{i,j} c_ij . d_ij
  /// summed over ALL ordered pairs (i,j) of occupied MOs (including
  /// i==j).
  static Eigen::MatrixXd RIKGradient(const Eigen::MatrixXd& occ_mo_coeffs,
                                      const AOBasis& auxbasis,
                                      const AOBasis& dftbasis);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTGRADIENT_H
