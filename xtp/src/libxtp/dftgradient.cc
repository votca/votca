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
#include "votca/xtp/screening_kernels.h"

namespace votca {
namespace xtp {

// Defined in libint2_derivative_calls.cc (no header).
using AOMatrixDerivative = std::array<Eigen::MatrixXd, 3>;
using ThreeCenterDerivative = std::array<std::vector<Eigen::MatrixXd>, 3>;
std::vector<AOMatrixDerivative> ComputeCoulombMetricDerivatives(
    const AOBasis& aobasis);
std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivatives(
    const AOBasis& auxbasis, const AOBasis& dftbasis);
// Density-contracted version of ComputeThreeCenterDerivatives, which stores
// the full tensor and is only usable for small systems.
std::vector<Eigen::MatrixXd> ComputeThreeCenterDerivativeContraction(
    const AOBasis& auxbasis, const AOBasis& dftbasis,
    const Eigen::MatrixXd& density);
// Occupied-MO-transformed three-centre derivatives, used by RIKGradient.
std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivativesMOTransformed(
    const AOBasis& auxbasis, const AOBasis& dftbasis,
    const Eigen::MatrixXd& occ_mo_coeffs);
// Per-atom variant, superseded by ComputeThreeCenterDerivativesMOTransformed.
ThreeCenterDerivative ComputeThreeCenterDerivativesForAtom(
    const AOBasis& auxbasis, const AOBasis& dftbasis, Index target_atom);
std::vector<Eigen::MatrixXd> ComputeThreeCenterIntegrals(
    const AOBasis& auxbasis, const AOBasis& dftbasis);

Eigen::MatrixXd DFTGradient::NuclearRepulsionDerivative(const QMMolecule& mol) {
  Index natoms = mol.size();
  Eigen::MatrixXd deriv = Eigen::MatrixXd::Zero(natoms, 3);

  // dE_nn/dR_A = sum_{B != A} Z_A Z_B * d(1/R_AB)/dR_A
  //            = -sum_{B != A} Z_A Z_B (R_A - R_B) / |R_A - R_B|^3
  //
  // Returns the gradient dE/dR, not the force.
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

  // V_PQ = (P|Q), c = V^-1 d. Plain solve, without the eigenvalue
  // truncation of Pseudo_InvSqrt used in the SCF; may need revisiting for
  // nearly linearly dependent aux bases.
  AOCoulomb aocoulomb;
  aocoulomb.Fill(auxbasis);
  const Eigen::MatrixXd& V = aocoulomb.Matrix();
  // V is positive definite; Cholesky, LDLT as before if it is not
  const CholeskyFactor V_chol(V);
  Eigen::VectorXd c = V_chol.ok() ? Eigen::VectorXd(V_chol.Solve(d))
                                  : Eigen::VectorXd(V.ldlt().solve(d));

  // ddP_dR[a](xyz, p) = sum_{mu,nu} P_{mu,nu} d(mu,nu|p)/dR_a,xyz, contracted
  // on the fly (the full derivative tensor does not fit in memory for
  // realistic systems).
  std::vector<Eigen::MatrixXd> ddP_dR =
      ComputeThreeCenterDerivativeContraction(auxbasis, dftbasis, density);
  std::vector<AOMatrixDerivative> dV =
      ComputeCoulombMetricDerivatives(auxbasis);

  // dE_J/dR = sum_P c_P d(d_P)/dR - 1/2 sum_PQ c_P c_Q d(V_PQ)/dR with
  // d(d_P)/dR = sum_{mu,nu} P_{mu,nu} d(mu,nu|P)/dR at fixed density (see
  // dftgradient.h).
  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      double term1 = ddP_dR[a].row(xyz).dot(c);
      double term2 = 0.5 * c.dot(dV[a][xyz] * c);
      grad(a, xyz) = term1 - term2;
    }
  }
  return grad;
}

Eigen::MatrixXd DFTGradient::RIKGradient(const Eigen::MatrixXd& occ_mo_coeffs,
                                         const AOBasis& auxbasis,
                                         const AOBasis& dftbasis) {
  Index natoms = static_cast<Index>(dftbasis.getFuncPerAtom().size());
  Index n_aux_bf = auxbasis.AOBasisSize();
  Index nocc = occ_mo_coeffs.cols();

  std::vector<Eigen::MatrixXd> tensor =
      ComputeThreeCenterIntegrals(auxbasis, dftbasis);

  AOCoulomb aocoulomb;
  aocoulomb.Fill(auxbasis);
  const Eigen::MatrixXd& V = aocoulomb.Matrix();
  // V is positive definite; Cholesky, LDLT as before if it is not
  const CholeskyFactor V_chol(V);
  Eigen::LDLT<Eigen::MatrixXd> V_ldlt;
  if (!V_chol.ok()) {
    V_ldlt.compute(V);
  }

  std::vector<AOMatrixDerivative> dV =
      ComputeCoulombMetricDerivatives(auxbasis);

  // d3c_mo[a][xyz][p](i, j) = C_i^T d(tensor[p])/dR_a,xyz C_j, all atoms in
  // one pass; memory natoms * 3 * Naux * nocc^2.
  std::vector<ThreeCenterDerivative> d3c_mo =
      ComputeThreeCenterDerivativesMOTransformed(auxbasis, dftbasis,
                                                 occ_mo_coeffs);

  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);

  // d_ij(P) = C_i^T tensor[P] C_j, c_ij = V^-1 d_ij over all ordered pairs
  // of occupied MOs i, j. The exchange energy is
  //   E_K = -2 * sum_{i,j} [0.5 * c_ij . d_ij] = -sum_{i,j} c_ij . d_ij,
  // equal to that of ERIs::CalculateEXX_mos.
  //
  // tensor[p] * C is formed once per p, so each d_ij(p) is an O(nao) dot
  // product.
  std::vector<Eigen::MatrixXd> tensor_half(n_aux_bf);
  for (Index p = 0; p < n_aux_bf; ++p) {
    tensor_half[p] = tensor[p] * occ_mo_coeffs;  // (nao, nocc)
  }
  std::vector<std::vector<Eigen::VectorXd>> c_ij(
      nocc, std::vector<Eigen::VectorXd>(nocc));
  for (Index i = 0; i < nocc; ++i) {
    for (Index j = 0; j < nocc; ++j) {
      Eigen::VectorXd d(n_aux_bf);
      for (Index p = 0; p < n_aux_bf; ++p) {
        d(p) = occ_mo_coeffs.col(i).dot(tensor_half[p].col(j));
      }
      c_ij[i][j] = V_chol.ok() ? Eigen::VectorXd(V_chol.Solve(d))
                               : Eigen::VectorXd(V_ldlt.solve(d));
    }
  }

  // CalculateEXX_mos uses the symmetric V^-1/2 fit; the V^-1 fit used here
  // gives the same energy (|V^-1/2 x|^2 = x^T V^-1 x for positive-definite
  // V), so no derivative of V^-1/2 is needed. Both MO indices must be
  // transformed: a half-transformed (one MO, one AO index) expression does
  // not reproduce the energy.

#pragma omp parallel for
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      double energy_term = 0.0;
      double metric_term = 0.0;
      for (Index i = 0; i < nocc; ++i) {
        for (Index j = 0; j < nocc; ++j) {
          const Eigen::VectorXd& c = c_ij[i][j];
          // d(d_ij(p))/dR_a,xyz
          Eigen::VectorXd dd(n_aux_bf);
          for (Index p = 0; p < n_aux_bf; ++p) {
            dd(p) = d3c_mo[a][xyz][p](i, j);
          }
          energy_term += c.dot(dd);
          metric_term += c.dot(dV[a][xyz] * c);
        }
      }
      // dE_K/dR = -2 sum_ij [c_ij . d(d_ij)/dR - 0.5 c_ij . dV/dR c_ij]
      grad(a, xyz) = -2.0 * (energy_term - 0.5 * metric_term);
    }
  }
  return grad;
}

}  // namespace xtp
}  // namespace votca
