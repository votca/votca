/*
 *            Copyright 2009-2021 The VOTCA Development Team
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

 #include "votca/xtp/rpa_uks.h"

#include "votca/xtp/aomatrix.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"

namespace votca {
namespace xtp {

void RPA_UKS::UpdateRPAInputEnergies(const Eigen::VectorXd& dftenergies_alpha,
                                     const Eigen::VectorXd& dftenergies_beta,
                                     const Eigen::VectorXd& gwaenergies_alpha,
                                     const Eigen::VectorXd& gwaenergies_beta,
                                     Index qpmin) {
  // Total number of orbitals retained in the RPA window.
  const Index rpatotal = rpamax_ - rpamin_ + 1;

  // Start from the DFT energies in the configured RPA window for each spin.
  energies_alpha_ = dftenergies_alpha.segment(rpamin_, rpatotal);
  energies_beta_ = dftenergies_beta.segment(rpamin_, rpatotal);

  // Number of GW-corrected states available in the qp window for each spin.
  const Index gwsize_alpha = Index(gwaenergies_alpha.size());
  const Index gwsize_beta = Index(gwaenergies_beta.size());

  // Replace DFT energies by GW energies where available.
  energies_alpha_.segment(qpmin - rpamin_, gwsize_alpha) = gwaenergies_alpha;
  energies_beta_.segment(qpmin - rpamin_, gwsize_beta) = gwaenergies_beta;

  // Outside the explicitly corrected qp window, shift remaining occupied and
  // virtual states by the largest observed correction in the corresponding
  // sector. This follows the same idea as the restricted implementation, but
  // is applied independently to alpha and beta channels.
  ShiftUncorrectedEnergies(energies_alpha_, dftenergies_alpha, homo_alpha_,
                           qpmin, gwsize_alpha);
  ShiftUncorrectedEnergies(energies_beta_, dftenergies_beta, homo_beta_, qpmin,
                           gwsize_beta);
}

void RPA_UKS::ShiftUncorrectedEnergies(Eigen::VectorXd& energies,
                                       const Eigen::VectorXd& dftenergies,
                                       Index homo, Index qpmin,
                                       Index gwsize) {
  const Index lumo = homo + 1;
  const Index qpmax = qpmin + gwsize - 1;

  // Largest absolute correction seen for occupied states inside the explicitly
  // corrected GW window.
  const double max_correction_occ =
      getMaxCorrection(energies, dftenergies, qpmin, homo);

  // Largest absolute correction seen for virtual states inside the explicitly
  // corrected GW window.
  const double max_correction_virt =
      getMaxCorrection(energies, dftenergies, lumo, qpmax);

  // The local "energies" vector is indexed relative to rpamin_, not relative to
  // the absolute orbital numbering. Therefore qpmin and qpmax must first be
  // shifted by rpamin_ when selecting head/tail segments.
  energies.head(qpmin - rpamin_).array() -= max_correction_occ;
  energies.tail(rpamax_ - qpmax).array() += max_correction_virt;
}

double RPA_UKS::getMaxCorrection(const Eigen::VectorXd& energies,
                                 const Eigen::VectorXd& dftenergies, Index min,
                                 Index max) const {
  if (max < min) {
    return 0.0;
  }

  const Index range = max - min + 1;

  // Compare the current working energies (possibly partly GW-corrected) to the
  // original DFT energies over the requested orbital range and return the
  // largest absolute deviation.
  const Eigen::VectorXd corrections =
      energies.segment(min - rpamin_, range) - dftenergies.segment(min, range);
  return corrections.cwiseAbs().maxCoeff();
}

template <bool imag>
Eigen::MatrixXd RPA_UKS::calculate_epsilon(double frequency) const {
  // The dielectric matrix lives in the auxiliary basis.
  const Index size = Mmn_.alpha.auxsize();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(size, size);

  // Accumulate the independent-particle polarizability contribution for one
  // spin channel. The final dielectric matrix is the sum of alpha and beta
  // contributions.
  auto accumulate_spin = [&](const TCMatrix_gwbse& Mmn,
                             const Eigen::VectorXd& energies, Index homo) {
    const Index lumo = homo + 1;

    // Number of occupied and unoccupied orbitals in the selected RPA window for
    // this spin channel.
    const Index n_occ = lumo - rpamin_;
    const Index n_unocc = rpamax_ - lumo + 1;

    if (n_occ <= 0 || n_unocc <= 0) {
      return;
    }

    const double freq2 = frequency * frequency;
    const double eta2 = eta_ * eta_;

    OpenMP_CUDA transform;
    transform.createTemporaries(n_unocc, size);

#pragma omp parallel
    {
      const Index threadid = OPENMP::getThreadId();

#pragma omp for schedule(dynamic)
      for (Index m_level = 0; m_level < n_occ; m_level++) {
        // m_level labels an occupied state within the local [rpamin_, rpamax_]
        // window. The associated virtual states are taken from the bottom part
        // of the same local energy ladder.
        const double qp_energy_m = energies(m_level);

        // Mmn[m_level] stores the three-center objects for a fixed occupied
        // state m against all n in the selected orbital window. bottomRows(...)
        // extracts the virtual block n = lumo ... rpamax for this spin.
        Eigen::MatrixXd Mmn_RPA = Mmn[m_level].bottomRows(n_unocc);
        transform.PushMatrix(Mmn_RPA, threadid);

        // Bare particle-hole energy differences Delta_e = e_a - e_i for this
        // fixed occupied state i = m_level and all virtual states a.
        const Eigen::ArrayXd deltaE =
            energies.tail(n_unocc).array() - qp_energy_m;

        Eigen::VectorXd denom;

        if (imag) {
          // Imaginary-frequency expression:
          //
          //   contribution ~ 2 * Delta_e / (Delta_e^2 + w^2)
          //
          // This matches the algebraic structure used in the restricted code.
          // The crucial difference is that here alpha and beta channels are
          // added explicitly instead of being folded into a closed-shell
          // degeneracy factor.
          denom = 2.0 * deltaE / (deltaE.square() + freq2);
        } else {
          // Real-frequency expression with finite broadening eta:
          //
          //   2 * [ (Delta_e - w)/((Delta_e - w)^2 + eta^2)
          //       + (Delta_e + w)/((Delta_e + w)^2 + eta^2) ]
          //
          // Again, the explicit spin resolution is handled by summing over the
          // two channels separately.
          Eigen::ArrayXd deltEf = deltaE - frequency;
          Eigen::ArrayXd sum = deltEf / (deltEf.square() + eta2);
          deltEf = deltaE + frequency;
          sum += deltEf / (deltEf.square() + eta2);
          denom = sum;
        }

        // Contract the virtual-index weights back to the auxiliary basis. This
        // is the same optimized RI reduction pattern as in the restricted code.
        transform.A_TDA(denom, threadid);
      }
    }

    // Add this spin-channel contribution to the final dielectric matrix.
    result += transform.getReductionVar();
  };

  // Total independent-particle response is the sum of alpha and beta channels.
  accumulate_spin(Mmn_.alpha, energies_alpha_, homo_alpha_);
  accumulate_spin(Mmn_.beta, energies_beta_, homo_beta_);

  // epsilon = 1 - v chi0  in the present RI convention ultimately appears here
  // as an identity contribution on the auxiliary-space diagonal added to the
  // accumulated response matrix, consistent with the restricted implementation.
  result.diagonal().array() += 1.0;
  return result;
}

template Eigen::MatrixXd RPA_UKS::calculate_epsilon<true>(double frequency) const;
template Eigen::MatrixXd RPA_UKS::calculate_epsilon<false>(double frequency) const;

Eigen::MatrixXd RPA_UKS::calculate_epsilon_r(
    std::complex<double> frequency) const {
  const Index size = Mmn_.alpha.auxsize();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(size, size);

  // Same structure as above, but for the real part evaluated at complex
  // frequency z = w + i*gamma.
  auto accumulate_spin = [&](const TCMatrix_gwbse& Mmn,
                             const Eigen::VectorXd& energies, Index homo) {
    const Index lumo = homo + 1;
    const Index n_occ = lumo - rpamin_;
    const Index n_unocc = rpamax_ - lumo + 1;

    if (n_occ <= 0 || n_unocc <= 0) {
      return;
    }

    OpenMP_CUDA transform;
    transform.createTemporaries(n_unocc, size);

#pragma omp parallel
    {
      const Index threadid = OPENMP::getThreadId();

#pragma omp for schedule(dynamic)
      for (Index m_level = 0; m_level < n_occ; m_level++) {
        const double qp_energy_m = energies(m_level);

        Eigen::MatrixXd Mmn_RPA = Mmn[m_level].bottomRows(n_unocc);
        transform.PushMatrix(Mmn_RPA, threadid);

        const Eigen::ArrayXd deltaE =
            energies.tail(n_unocc).array() - qp_energy_m;

        // Terms corresponding to Re[ 1 / (w + i*gamma - Delta_e) ] and its
        // partner at -Delta_e, written in the same real-valued form as in the
        // restricted implementation.
        const Eigen::ArrayXd deltaEm = frequency.real() - deltaE;
        const Eigen::ArrayXd deltaEp = frequency.real() + deltaE;

        const double sigma_1 = std::pow(frequency.imag() + eta_, 2);
        const double sigma_2 = std::pow(frequency.imag() - eta_, 2);

        const Eigen::VectorXd chi =
            deltaEm * (deltaEm.cwiseAbs2() + sigma_1).cwiseInverse() -
            deltaEp * (deltaEp.cwiseAbs2() + sigma_2).cwiseInverse();

        transform.A_TDA(chi, threadid);
      }
    }

    // The overall prefactor follows the existing restricted convention. The
    // spin resolution still enters explicitly through separate alpha/beta sums.
    result -=  transform.getReductionVar();
  };

  accumulate_spin(Mmn_.alpha, energies_alpha_, homo_alpha_);
  accumulate_spin(Mmn_.beta, energies_beta_, homo_beta_);

  result.diagonal().array() += 1.0;
  return result;
}

RPA_UKS::rpa_eigensolution RPA_UKS::Diagonalize_H2p() const {
  // AmB contains the bare particle-hole energy differences and is diagonal in
  // the particle-hole basis.
  Eigen::VectorXd AmB = Calculate_H2p_AmB();

  // ApB adds the Coulomb coupling between particle-hole states, including
  // alpha-alpha, beta-beta, and mixed alpha-beta / beta-alpha blocks.
  Eigen::MatrixXd ApB = Calculate_H2p_ApB();

  RPA_UKS::rpa_eigensolution sol;
  sol.ERPA_correlation = -0.25 * (ApB.trace() + AmB.sum());

  // Build the symmetrized Hermitian matrix
  //
  //   C = sqrt(AmB) (ApB) sqrt(AmB)
  //
  // exactly as done in the restricted code path.
  Eigen::MatrixXd& C = ApB;
  C.applyOnTheLeft(AmB.cwiseSqrt().asDiagonal());
  C.applyOnTheRight(AmB.cwiseSqrt().asDiagonal());

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es = Diagonalize_H2p_C(C);

  sol.omega = Eigen::VectorXd::Zero(es.eigenvalues().size());
  sol.omega = es.eigenvalues().cwiseSqrt();
  sol.ERPA_correlation += 0.5 * sol.omega.sum();

  XTP_LOG(Log::info, log_) << TimeStamp()
                           << " Lowest neutral excitation energy (eV): "
                           << tools::conv::hrt2ev * sol.omega.minCoeff()
                           << std::flush;

  XTP_LOG(Log::error, log_)
      << TimeStamp()
      << " RPA correlation energy (Hartree): " << sol.ERPA_correlation
      << std::flush;

  const Index rpasize = Index(AmB.size());
  sol.XpY = Eigen::MatrixXd(rpasize, rpasize);

  // Reconstruct X+Y eigenvectors from the symmetrized eigenvectors.
  const Eigen::VectorXd AmB_sqrt = AmB.cwiseSqrt();
  const Eigen::VectorXd Omega_sqrt_inv = sol.omega.cwiseSqrt().cwiseInverse();
  for (Index s = 0; s < rpasize; s++) {
    sol.XpY.col(s) =
        Omega_sqrt_inv(s) * AmB_sqrt.cwiseProduct(es.eigenvectors().col(s));
  }

  return sol;
}

Eigen::VectorXd RPA_UKS::Calculate_H2p_AmB() const {
  const Index lumo_alpha = homo_alpha_ + 1;
  const Index lumo_beta = homo_beta_ + 1;

  const Index n_occ_alpha = lumo_alpha - rpamin_;
  const Index n_occ_beta = lumo_beta - rpamin_;
  const Index n_unocc_alpha = rpamax_ - lumo_alpha + 1;
  const Index n_unocc_beta = rpamax_ - lumo_beta + 1;

  // Total number of alpha and beta particle-hole excitations in the selected
  // window.
  const Index size_alpha = n_occ_alpha * n_unocc_alpha;
  const Index size_beta = n_occ_beta * n_unocc_beta;

  Eigen::VectorXd AmB = Eigen::VectorXd::Zero(size_alpha + size_beta);

  // alpha block:
  // excitation basis index enumerates all alpha (v -> c) combinations
  vc2index vc_alpha(0, 0, n_unocc_alpha);
  for (Index v = 0; v < n_occ_alpha; v++) {
    const Index i = vc_alpha.I(v, 0);
    AmB.segment(i, n_unocc_alpha) =
        energies_alpha_.segment(n_occ_alpha, n_unocc_alpha).array() -
        energies_alpha_(v);
  }

  // beta block:
  // appended after all alpha excitations
  vc2index vc_beta(0, 0, n_unocc_beta);
  for (Index v = 0; v < n_occ_beta; v++) {
    const Index i = size_alpha + vc_beta.I(v, 0);
    AmB.segment(i, n_unocc_beta) =
        energies_beta_.segment(n_occ_beta, n_unocc_beta).array() -
        energies_beta_(v);
  }

  return AmB;
}

Eigen::MatrixXd RPA_UKS::Calculate_H2p_ApB() const {
  const Index lumo_alpha = homo_alpha_ + 1;
  const Index lumo_beta = homo_beta_ + 1;

  const Index n_occ_alpha = lumo_alpha - rpamin_;
  const Index n_occ_beta = lumo_beta - rpamin_;
  const Index n_unocc_alpha = rpamax_ - lumo_alpha + 1;
  const Index n_unocc_beta = rpamax_ - lumo_beta + 1;

  const Index size_alpha = n_occ_alpha * n_unocc_alpha;
  const Index size_beta = n_occ_beta * n_unocc_beta;

  Eigen::MatrixXd ApB =
      Eigen::MatrixXd::Zero(size_alpha + size_beta, size_alpha + size_beta);

  vc2index vc_alpha(0, 0, n_unocc_alpha);
  vc2index vc_beta(0, 0, n_unocc_beta);

  // alpha-alpha block:
  // Coulomb coupling between alpha particle-hole excitations.
#pragma omp parallel for schedule(guided)
  for (Index v2 = 0; v2 < n_occ_alpha; v2++) {
    const Index i2 = vc_alpha.I(v2, 0);
    const Eigen::MatrixXd Mmn_v2T =
        Mmn_.alpha[v2].middleRows(n_occ_alpha, n_unocc_alpha).transpose();

    for (Index v1 = v2; v1 < n_occ_alpha; v1++) {
      const Index i1 = vc_alpha.I(v1, 0);

      // Only a single factor 2.0 is kept here from the algebraic A+B
      // structure. The closed-shell spin degeneracy factor of the restricted
      // implementation is not present anymore; spin is represented explicitly
      // by separate alpha and beta blocks.
      ApB.block(i1, i2, n_unocc_alpha, n_unocc_alpha) =
          1.0 * Mmn_.alpha[v1].middleRows(n_occ_alpha, n_unocc_alpha) * Mmn_v2T;
    }
  }

  // beta-beta block:
#pragma omp parallel for schedule(guided)
  for (Index v2 = 0; v2 < n_occ_beta; v2++) {
    const Index i2 = size_alpha + vc_beta.I(v2, 0);
    const Eigen::MatrixXd Mmn_v2T =
        Mmn_.beta[v2].middleRows(n_occ_beta, n_unocc_beta).transpose();

    for (Index v1 = v2; v1 < n_occ_beta; v1++) {
      const Index i1 = size_alpha + vc_beta.I(v1, 0);
      ApB.block(i1, i2, n_unocc_beta, n_unocc_beta) =
          1.0 * Mmn_.beta[v1].middleRows(n_occ_beta, n_unocc_beta) * Mmn_v2T;
    }
  }

  // alpha-beta block:
  // In RPA, Coulomb coupling acts on total density fluctuations and therefore
  // couples alpha and beta particle-hole sectors as well.
#pragma omp parallel for schedule(guided)
  for (Index v_beta = 0; v_beta < n_occ_beta; v_beta++) {
    const Index i_beta = size_alpha + vc_beta.I(v_beta, 0);
    const Eigen::MatrixXd Mmn_beta_T =
        Mmn_.beta[v_beta].middleRows(n_occ_beta, n_unocc_beta).transpose();

    for (Index v_alpha = 0; v_alpha < n_occ_alpha; v_alpha++) {
      const Index i_alpha = vc_alpha.I(v_alpha, 0);
      ApB.block(i_alpha, i_beta, n_unocc_alpha, n_unocc_beta) =
          1.0 *
          Mmn_.alpha[v_alpha].middleRows(n_occ_alpha, n_unocc_alpha) *
          Mmn_beta_T;
    }
  }

  // Fill the opposite mixed block and the missing lower-triangular parts by
  // Hermitian symmetry.
  ApB.template triangularView<Eigen::StrictlyUpper>() =
      ApB.transpose().template triangularView<Eigen::StrictlyUpper>();

  // Add the diagonal bare transition energies.
  ApB.diagonal() += Calculate_H2p_AmB();

  return ApB;
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> RPA_UKS::Diagonalize_H2p_C(
    const Eigen::MatrixXd& C) const {
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Diagonalizing two-particle Hamiltonian "
      << std::flush;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Diagonalization done " << std::flush;

  const double minCoeff = es.eigenvalues().minCoeff();
  if (minCoeff <= 0.0) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Detected non-positive eigenvalue: " << minCoeff
        << std::flush;
    throw std::runtime_error("Detected non-positive eigenvalue.");
  }

  return es;
}

}  // namespace xtp
}  // namespace votca