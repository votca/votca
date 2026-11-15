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

// Contracted RI integrals and integral derivatives for analytic gradients
// (defined in libint2_derivative_calls.cc). Nothing of size
// Naux * Nao^2 or natoms * Naux^2 is stored: three-centre integrals and
// their derivatives are contracted per aux shell as they are computed, the
// metric derivatives per shell pair.
//
// Any RI gradient term of the form
//   sum_P sum_{mu nu} Gamma^P_{mu nu} d(mu nu|P)/dR  -  1/2 sum_PQ G_PQ
//   d(P|Q)/dR
// (RI-J, RI-K, and the screened-exchange terms of BSE gradients) is
// assembled from ContractThreeCenterDerivatives and
// ContractCoulombMetricDerivatives.

#pragma once
#ifndef VOTCA_XTP_RIDERIVATIVES_H
#define VOTCA_XTP_RIDERIVATIVES_H

// Standard includes
#include <functional>

// Local VOTCA includes
#include "aobasis.h"
#include "eigen.h"

namespace votca {
namespace xtp {

/// Called with an aux function index P and the (Nao x Nao) matrix
/// (mu nu|P). Called concurrently from several threads, for different P.
using ThreeCenterVisitor = std::function<void(Index, const Eigen::MatrixXd&)>;

/// Fills Gamma (resized by the callee to Nao x Nao) with the AO three-index
/// density Gamma^P of aux function P. Only the symmetric part enters, since
/// (mu nu|P) is symmetric in mu, nu. Called concurrently from several
/// threads, for different P.
using AuxDensityFn = std::function<void(Index, Eigen::MatrixXd&)>;

/// Calls use(P, (mu nu|P)) for every aux function P, in parallel over aux
/// shells; memory one aux shell of Nao x Nao matrices per thread.
void ForEachThreeCenterIntegral(const AOBasis& auxbasis,
                                const AOBasis& dftbasis,
                                const ThreeCenterVisitor& use);

/// grad(a, xyz) = sum_P sum_{mu nu} Gamma^P_{mu nu} d(mu nu|P)/dR_a,xyz,
/// (Natoms x 3). Gamma^P is requested once per aux function.
Eigen::MatrixXd ContractThreeCenterDerivatives(const AOBasis& auxbasis,
                                               const AOBasis& dftbasis,
                                               const AuxDensityFn& gamma);

/// grad(a, xyz) = sum_PQ G_PQ d(P|Q)/dR_a,xyz for a symmetric (Naux x Naux)
/// G, (Natoms x 3). Natoms is that of the aux basis.
Eigen::MatrixXd ContractCoulombMetricDerivatives(const AOBasis& auxbasis,
                                                 const Eigen::MatrixXd& G);

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_RIDERIVATIVES_H
