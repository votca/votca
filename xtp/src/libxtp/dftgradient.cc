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

// Local VOTCA includes
#include "votca/xtp/dftgradient.h"
#include "votca/xtp/aomatrix.h"

namespace votca {
namespace xtp {

// Defined in libint2_derivative_calls.cc, not yet in any header (same
// STATUS noted throughout that file and libint2_derivative_calls.cc
// itself) -- forward declared here rather than adding a new header at
// this stage, consistent with how these functions have been consumed
// from test files so far.
using AOMatrixDerivative = std::array<Eigen::MatrixXd, 3>;
using ThreeCenterDerivative = std::array<std::vector<Eigen::MatrixXd>, 3>;
std::vector<AOMatrixDerivative> ComputeCoulombMetricDerivatives(
    const AOBasis& aobasis);
std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivatives(
    const AOBasis& auxbasis, const AOBasis& dftbasis);
std::vector<Eigen::MatrixXd> ComputeThreeCenterIntegrals(
    const AOBasis& auxbasis, const AOBasis& dftbasis);

Eigen::MatrixXd DFTGradient::NuclearRepulsionDerivative(
    const QMMolecule& mol) {
  Index natoms = mol.size();
  Eigen::MatrixXd deriv = Eigen::MatrixXd::Zero(natoms, 3);

  // dE_nn/dR_A = sum_{B != A} Z_A Z_B * d(1/R_AB)/dR_A
  //            = -sum_{B != A} Z_A Z_B (R_A - R_B) / |R_A - R_B|^3
  //
  // NOTE: this returns the GRADIENT dE/dR, not the force -dE/dR -- the
  // sign convention for the final assembled total gradient (which this
  // feeds into) is a decision for whatever code combines this with the
  // RI-J and XC terms, not fixed here. Keep this consistent when
  // validating against finite differences: a finite difference of the
  // energy directly gives dE/dR, matching this function's output as-is,
  // with no extra sign flip needed.
  for (Index a = 0; a < natoms; ++a) {
    double Za = static_cast<double>(mol[a].getNuccharge());
    const Eigen::Vector3d& Ra = mol[a].getPos();
    Eigen::Vector3d sum = Eigen::Vector3d::Zero();
    for (Index b = 0; b < natoms; ++b) {
      if (b == a) {
        continue;
      }
      double Zb = static_cast<double>(mol[b].getNuccharge());
      const Eigen::Vector3d& Rb = mol[b].getPos();
      Eigen::Vector3d Rab_vec = Ra - Rb;
      double Rab = Rab_vec.norm();
      sum += Za * Zb * Rab_vec / (Rab * Rab * Rab);
    }
    deriv.row(a) = -sum.transpose();
  }
  return deriv;
}

Eigen::MatrixXd DFTGradient::RIJGradient(const Eigen::MatrixXd& density,
                                          const AOBasis& auxbasis,
                                          const AOBasis& dftbasis) {
  Index natoms = static_cast<Index>(dftbasis.getFuncPerAtom().size());
  Index n_aux_bf = auxbasis.AOBasisSize();

  // d_P = sum_{mu,nu} P_{mu,nu} (mu,nu|P)
  std::vector<Eigen::MatrixXd> tensor =
      ComputeThreeCenterIntegrals(auxbasis, dftbasis);
  Eigen::VectorXd d(n_aux_bf);
  for (Index p = 0; p < n_aux_bf; ++p) {
    d(p) = (density.array() * tensor[p].array()).sum();
  }

  // V_PQ = (P|Q), c = V^-1 d. Using a plain SPD solve here rather than
  // the eigenvalue-truncated Pseudo_InvSqrt approach AOCoulomb also
  // offers (used in production for numerical stability against
  // near-linear-dependence in the aux basis) -- adequate for this
  // validation-scale use, but worth revisiting if this is ever used on
  // a genuinely large or near-linearly-dependent aux basis.
  AOCoulomb aocoulomb;
  aocoulomb.Fill(auxbasis);
  const Eigen::MatrixXd& V = aocoulomb.Matrix();
  Eigen::VectorXd c = V.ldlt().solve(d);

  // Derivative tensors -- both already validated (finite-difference
  // tested) in test_aoderivatives.cc.
  std::vector<ThreeCenterDerivative> d3c =
      ComputeThreeCenterDerivatives(auxbasis, dftbasis);
  std::vector<AOMatrixDerivative> dV = ComputeCoulombMetricDerivatives(auxbasis);

  // dE_J/dR = sum_P c_P d(d_P)/dR - 1/2 sum_PQ c_P c_Q d(V_PQ)/dR
  // d(d_P)/dR = sum_{mu,nu} P_{mu,nu} d(mu,nu|P)/dR (density held fixed --
  // see the "IMPORTANT" note on this function in dftgradient.h for why
  // that's valid regardless of whether density is a converged SCF
  // density or an arbitrary fixed matrix).
  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      double term1 = 0.0;
      for (Index p = 0; p < n_aux_bf; ++p) {
        term1 += c(p) * (density.array() * d3c[a][xyz][p].array()).sum();
      }
      double term2 = 0.5 * c.dot(dV[a][xyz] * c);
      grad(a, xyz) = term1 - term2;
    }
  }
  return grad;
}

}  // namespace xtp
}  // namespace votca
