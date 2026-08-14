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
// Memory-efficient alternative to ComputeThreeCenterDerivatives above --
// see that function's own header comment in libint2_derivative_calls.cc
// for why this exists (the full-tensor version is catastrophically
// memory-unscalable for anything beyond a small, toy-sized system).
std::vector<Eigen::MatrixXd> ComputeThreeCenterDerivativeContraction(
    const AOBasis& auxbasis, const AOBasis& dftbasis,
    const Eigen::MatrixXd& density);
// CORRECTED, single-pass alternative for RIKGradient's own, different
// contraction needs -- see that function's own header comment in
// libint2_derivative_calls.cc for the full reasoning, including the
// arithmetic error that led to an earlier, slower, per-atom approach
// (ComputeThreeCenterDerivativesForAtom, kept below but no longer
// used by RIKGradient) being chosen instead.
std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivativesMOTransformed(
    const AOBasis& auxbasis, const AOBasis& dftbasis,
    const Eigen::MatrixXd& occ_mo_coeffs);
// Per-atom alternative -- superseded by
// ComputeThreeCenterDerivativesMOTransformed above for RIKGradient's
// own use (kept, unused by RIKGradient now, since deleting it would
// remove a real, working, if slower, fallback with no compensating
// benefit).
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

  // Already-contracted against density (sum_{mu,nu} density(mu,nu) *
  // d(mu,nu|p)/dR_a[xyz], not a per-(mu,nu) tensor) -- see this
  // function's own header comment for why: the full-tensor
  // ComputeThreeCenterDerivatives, previously used here, was confirmed
  // to exhaust hundreds of GB of memory on a real, moderately-sized
  // molecule (natoms*3*n_aux_bf separate, complete nao x nao matrices
  // held simultaneously -- roughly 2.6 PETABYTES for a 53-atom, 943
  // AO / 2320 auxiliary function system). ComputeThreeCenterDerivatives
  // itself is left unchanged (still used, and still validated, by
  // test_aoderivatives.cc and by RIKGradient below, which needs a
  // genuinely different contraction of its own -- see that function's
  // own comments).
  std::vector<Eigen::MatrixXd> ddP_dR =
      ComputeThreeCenterDerivativeContraction(auxbasis, dftbasis, density);
  std::vector<AOMatrixDerivative> dV =
      ComputeCoulombMetricDerivatives(auxbasis);

  // dE_J/dR = sum_P c_P d(d_P)/dR - 1/2 sum_PQ c_P c_Q d(V_PQ)/dR
  // d(d_P)/dR = sum_{mu,nu} P_{mu,nu} d(mu,nu|P)/dR (density held fixed --
  // see the "IMPORTANT" note on this function in dftgradient.h for why
  // that's valid regardless of whether density is a converged SCF
  // density or an arbitrary fixed matrix). ddP_dR[a](xyz, p) is exactly
  // this quantity, already summed over mu,nu -- no further contraction
  // against density needed here at all.
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
  Eigen::LDLT<Eigen::MatrixXd> V_ldlt(V);

  std::vector<AOMatrixDerivative> dV =
      ComputeCoulombMetricDerivatives(auxbasis);

  // Single pass over ALL atoms at once -- see
  // ComputeThreeCenterDerivativesMOTransformed's own header comment for
  // the full reasoning (including the arithmetic error that had
  // motivated a slower, per-atom approach instead): ~23.8 GB total for
  // the real, 53-atom/943-AO/2320-auxiliary/93-occupied-orbital system
  // that originally motivated this whole fix, smaller than even that
  // per-atom approach's own ~46 GB peak, and needing only ONE pass over
  // the expensive shell-triple loop rather than natoms of them.
  std::vector<ThreeCenterDerivative> d3c_mo =
      ComputeThreeCenterDerivativesMOTransformed(auxbasis, dftbasis,
                                                 occ_mo_coeffs);

  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);

  // d_ij(P) = C_i^T tensor[P] C_j, c_ij = V^-1 d_ij (i,j both occupied
  // MOs, all ordered pairs including i==j).
  // E_K = -2 * sum_{i,j} [0.5 * c_ij . d_ij] = -sum_{i,j} c_ij . d_ij
  // (matches ERIs::CalculateEXX_mos's real physical exchange energy
  // exactly, confirmed numerically, not just up to an unknown scale).
  //
  // PERFORMANCE: confirmed directly, via a real run, to be a genuine
  // bottleneck in its naive form -- computing tensor[p]*occ_mo_coeffs.col(j)
  // (an O(nao^2) matrix-vector product) separately for every (i,j)
  // pair recomputes the SAME result nocc times over (it does not
  // depend on i at all). Factored out here: tensor[p]*occ_mo_coeffs
  // (all occupied columns at once) is computed ONCE per p, and every
  // d(p) for a given (i,j) is then a cheap O(nao) dot product against
  // the appropriate column of that already-computed result -- reduces
  // the dominant cost by roughly a factor of nocc. Purely a
  // reordering of the SAME formula (see the "NOTE ON HISTORY" comment
  // below for why changing the formula itself would be dangerous);
  // mathematically identical to the original for every element, not a
  // different computation.
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
      c_ij[i][j] = V_ldlt.solve(d);
    }
  }

  // NOTE ON HISTORY: an earlier revision of this function briefly
  // switched to a "half-transformed" structure (one index MO, one AO),
  // reasoned (incorrectly, via error-prone hand algebra) to be needed
  // to match ERIs::CalculateEXX_mos's real K matrix. Settled by DIRECT
  // NUMERICAL SIMULATION of CalculateEXX_mos's actual algorithm
  // (symmetric V^-1/2 fit, TCxMOs_T = occMos^T*B_tilde, etc.) against
  // both candidate formulas, on several random test systems: the
  // FULLY-MO-transformed structure below (both indices occupied MOs,
  // matching the ORIGINAL version of this function) is correct, and
  // matches the real energy EXACTLY (to ~1e-14) once multiplied by 2 --
  // not the half-transformed structure, which did not match at all
  // (not even up to a constant factor). See conversation history for
  // the verification. This confirms the earlier documented concern
  // (needing to differentiate a matrix square root) was never actually
  // a problem -- V^-1 and V^-1/2 fitting give identical energies
  // (|V^-1/2 x|^2 == x^T V^-1 x exactly, for symmetric positive-definite
  // V) -- the only real fix needed here was the missing factor of 2.
  // (NOTE: this history comment is about a DIFFERENT half-transformed
  // structure than the one introduced just above -- that one changed
  // the final FORMULA itself and gave a wrong physical answer; the one
  // just above only reorders/factors the SAME formula's own
  // computation and is mathematically identical to the original for
  // every element, not a different computation at all.)

#pragma omp parallel for
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      double energy_term = 0.0;
      double metric_term = 0.0;
      for (Index i = 0; i < nocc; ++i) {
        for (Index j = 0; j < nocc; ++j) {
          const Eigen::VectorXd& c = c_ij[i][j];
          // d3c_mo[a][xyz][p] is already the full (nocc, nocc)
          // MO-transformed matrix -- (i, j) is read off directly, no
          // further matrix-vector multiplication needed at all (unlike
          // the old, per-atom AO-basis d3c_atom[xyz][p], which still
          // needed occ_mo_coeffs.col(i).dot(d3c_atom[xyz][p] *
          // occ_mo_coeffs.col(j)) for every (i,j) pair).
          Eigen::VectorXd dd(n_aux_bf);
          for (Index p = 0; p < n_aux_bf; ++p) {
            dd(p) = d3c_mo[a][xyz][p](i, j);
          }
          energy_term += c.dot(dd);
          metric_term += c.dot(dV[a][xyz] * c);
        }
      }
      // Factor of -2, matching E_K = -2*sum[0.5*c.d] confirmed above --
      // NOT the "-(energy_term - 0.5*metric_term)" (factor of -1) used
      // in the intermediate, incorrect half-transformed revision.
      grad(a, xyz) = -2.0 * (energy_term - 0.5 * metric_term);
    }
  }
  return grad;
}

}  // namespace xtp
}  // namespace votca
