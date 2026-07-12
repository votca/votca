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
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTGRADIENT_H
