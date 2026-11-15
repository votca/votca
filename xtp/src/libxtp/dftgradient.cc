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
#include "votca/xtp/riderivatives.h"
#include "votca/xtp/screening_kernels.h"

namespace votca {
namespace xtp {

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

namespace {
// V^-1 B for the Coulomb metric V = (P|Q), Cholesky with LDLT as fallback
// (plain solve, without the eigenvalue truncation of Pseudo_InvSqrt used in
// the SCF; may need revisiting for nearly linearly dependent aux bases).
Eigen::MatrixXd SolveMetric(const AOBasis& auxbasis, const Eigen::MatrixXd& B) {
  AOCoulomb aocoulomb;
  aocoulomb.Fill(auxbasis);
  const Eigen::MatrixXd& V = aocoulomb.Matrix();
  const CholeskyFactor V_chol(V);
  if (V_chol.ok()) {
    return V_chol.Solve(B);
  }
  return V.ldlt().solve(B);
}
}  // namespace

Eigen::MatrixXd DFTGradient::RIJGradient(const Eigen::MatrixXd& density,
                                         const AOBasis& auxbasis,
                                         const AOBasis& dftbasis) {
  Index n_aux_bf = auxbasis.AOBasisSize();

  // d_P = sum_{mu,nu} D_{mu,nu} (mu,nu|P), c = V^-1 d
  Eigen::VectorXd d(n_aux_bf);
  ForEachThreeCenterIntegral(auxbasis, dftbasis,
                             [&density, &d](Index P, const Eigen::MatrixXd& I) {
                               d(P) = (density.array() * I.array()).sum();
                             });
  const Eigen::VectorXd c = SolveMetric(auxbasis, d);

  // dE_J/dR = sum_P c_P d(d_P)/dR - 1/2 sum_PQ c_P c_Q d(V_PQ)/dR with
  // d(d_P)/dR = sum_{mu,nu} D_{mu,nu} d(mu,nu|P)/dR at fixed density (see
  // dftgradient.h), i.e. Gamma^P = c_P D and G = c c^T.
  Eigen::MatrixXd grad = ContractThreeCenterDerivatives(
      auxbasis, dftbasis, [&density, &c](Index P, Eigen::MatrixXd& Gamma) {
        Gamma = c(P) * density;
      });
  grad -= 0.5 * ContractCoulombMetricDerivatives(auxbasis, c * c.transpose());
  return grad;
}

Eigen::MatrixXd DFTGradient::RIKGradient(const Eigen::MatrixXd& occ_mo_coeffs,
                                         const AOBasis& auxbasis,
                                         const AOBasis& dftbasis) {
  // E_K = -sum_{i,j} c_ij . d_ij with d_ij(P) = C_i^T (mu nu|P) C_j,
  // c_ij = V^-1 d_ij over all ordered pairs of occupied MOs, equal to the
  // exchange energy of ERIs::CalculateEXX_mos (symmetric V^-1/2 fit; the
  // V^-1 fit used here gives the same energy, |V^-1/2 x|^2 = x^T V^-1 x, so
  // no derivative of V^-1/2 is needed).
  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(occ_mo_coeffs.cols());
  return -RIPairGradient(occ_mo_coeffs, ones, occ_mo_coeffs, auxbasis,
                         dftbasis);
}

Eigen::MatrixXd DFTGradient::RIPairGradient(const Eigen::MatrixXd& U,
                                            const Eigen::VectorXd& w,
                                            const Eigen::MatrixXd& C,
                                            const AOBasis& auxbasis,
                                            const AOBasis& dftbasis) {
  Index n_aux_bf = auxbasis.AOBasisSize();
  Index nk = U.cols();
  Index ni = C.cols();

  // d_ki(P) as the rows of dT (Naux x nk*ni, column k + nk*i), c_ki =
  // V^-1 d_ki for all pairs in one multi-RHS solve. Memory O(Naux nk ni).
  Eigen::MatrixXd dT(n_aux_bf, nk * ni);
  ForEachThreeCenterIntegral(
      auxbasis, dftbasis,
      [&U, &C, &dT, nk, ni](Index P, const Eigen::MatrixXd& I) {
        Eigen::MatrixXd dki = U.transpose() * (I * C);
        dT.row(P) = Eigen::Map<const Eigen::RowVectorXd>(dki.data(), nk * ni);
      });
  const Eigen::MatrixXd cT = SolveMetric(auxbasis, dT);  // Naux x nk*ni
  dT.resize(0, 0);

  // dE/dR = sum_ki w_k [2 c_ki . d(d_ki)/dR - c_ki . dV/dR c_ki]
  //   with Gamma^P = U diag(w) c^P C^T (c^P the nk x ni matrix of c_ki(P))
  //   and G = sum_ki w_k c_ki c_ki^T.
  Eigen::MatrixXd grad = ContractThreeCenterDerivatives(
      auxbasis, dftbasis,
      [&U, &w, &C, &cT, nk, ni](Index P, Eigen::MatrixXd& Gamma) {
        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<>> cP(
            cT.data() + P, nk, ni, Eigen::InnerStride<>(cT.rows()));
        Gamma = U * (w.asDiagonal() * Eigen::MatrixXd(cP) * C.transpose());
      });
  Eigen::MatrixXd cw = cT;
  for (Index i = 0; i < ni; ++i) {
    cw.middleCols(i * nk, nk) *= w.asDiagonal();
  }
  const Eigen::MatrixXd G = cw * cT.transpose();
  return 2.0 * grad - ContractCoulombMetricDerivatives(auxbasis, G);
}

}  // namespace xtp
}  // namespace votca
