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

#ifndef VOTCA_XTP_EXCITONFOLLOW_H
#define VOTCA_XTP_EXCITONFOLLOW_H

// Local VOTCA includes
#include "eigen.h"
#include "orbitals.h"
#include "qmstate.h"

namespace votca {
namespace xtp {

/// Root following between geometries: overlaps of the BSE state ref_state
/// of ref (an earlier geometry) with all BSE states of type of cur,
///   O_i = |<X_ref|X_i> - <Y_ref|Y_i>|
/// with the excitons in their AO representation C_v X^T C_c^T and the
/// mixed AO overlap <chi(R_ref)|chi(R_cur)> between the two geometries (so
/// the moving basis is accounted for). Equals 1 for the same state at the
/// same geometry. Both orbitals need the same basis set and atoms.
Eigen::VectorXd ExcitonOverlaps(const Orbitals& ref, const QMState& ref_state,
                                const Orbitals& cur, QMStateType type);

/// The mixed AO overlap S_mk = <chi_m(R_a)|chi_k(R_b)> between the DFT basis
/// functions of two geometries of the same molecule.
Eigen::MatrixXd MixedAOOverlap(const Orbitals& a, const Orbitals& b);

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EXCITONFOLLOW_H
