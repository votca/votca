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
#include "votca/xtp/screening_kernels.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "votca/xtp/aomatrix.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"

namespace {

constexpr double kUKSApBPrefactor = 2.0;

}  // namespace

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

  // Outside the explicitly corrected qp window, shift the remaining occupied
  // and virtual states (OutOfWindowShift), independently for alpha and beta.
  ShiftUncorrectedEnergies(energies_alpha_, dftenergies_alpha, homo_alpha_,
                           qpmin, gwsize_alpha);
  ShiftUncorrectedEnergies(energies_beta_, dftenergies_beta, homo_beta_, qpmin,
                           gwsize_beta);

  InvalidateH2pCache();
}

void RPA_UKS::InvalidateH2pCache() {
  h2p_cached_ = false;
  h2p_solution_cache_ = rpa_eigensolution{};

  screening_cached_ = false;
  screening_omegas_.resize(0);
  screening_modes_.clear();
}

void RPA_UKS::GetCachedScreeningModes(
    const Eigen::VectorXd*& omegas,
    const std::vector<Eigen::VectorXd>*& modes) const {
  if (!screening_cached_) {
    BuildCachedScreeningModes();
  }
  omegas = &screening_omegas_;
  modes = &screening_modes_;
}

void RPA_UKS::BuildCachedScreeningModes() const {
  const rpa_eigensolution& rpa_solution = Diagonalize_H2p();
  const Eigen::MatrixXd& XpY = rpa_solution.XpY;
  const Eigen::VectorXd& omegas = rpa_solution.omega;

  const Index lumo_alpha = homo_alpha_ + 1;
  const Index lumo_beta = homo_beta_ + 1;

  const Index n_occ_alpha = lumo_alpha - rpamin_;
  const Index n_occ_beta = lumo_beta - rpamin_;
  const Index n_unocc_alpha = rpamax_ - homo_alpha_;
  const Index n_unocc_beta = rpamax_ - homo_beta_;

  const Index size_alpha = n_occ_alpha * n_unocc_alpha;
  const Index auxsize = Mmn_.alpha.auxsize();
  const Index nmodes = Index(omegas.size());

  // modes Z = sum_v M_v,virt^T (X+Y)_v over both spins, as two GEMMs over
  // the stacked virtual rows of the hole slices
  auto stacked_rows = [&](const TCMatrix_gwbse& Mmn, Index n_occ,
                          Index n_unocc) {
    vc2index vc(0, 0, n_unocc);
    Eigen::MatrixXd M(n_occ * n_unocc, auxsize);
#pragma omp parallel for schedule(dynamic)
    for (Index v = 0; v < n_occ; v++) {
      M.middleRows(vc.I(v, 0), n_unocc) = Mmn[v].middleRows(n_occ, n_unocc);
    }
    return M;
  };
  const Index size_beta = n_occ_beta * n_unocc_beta;
  Eigen::MatrixXd Z =
      stacked_rows(Mmn_.alpha, n_occ_alpha, n_unocc_alpha).transpose() *
      XpY.topRows(size_alpha);
  Z.noalias() += stacked_rows(Mmn_.beta, n_occ_beta, n_unocc_beta).transpose() *
                 XpY.middleRows(size_alpha, size_beta);

  std::vector<Eigen::VectorXd> all_modes;
  all_modes.reserve(nmodes);
  std::vector<double> all_norms;
  all_norms.reserve(nmodes);
  double max_norm = 0.0;
  for (Index s = 0; s < nmodes; s++) {
    const double norm = Z.col(s).norm();
    max_norm = std::max(max_norm, norm);
    all_modes.push_back(Z.col(s));
    all_norms.push_back(norm);
  }

  const double tol = 1e-10 * std::max(1.0, max_norm);

  std::vector<Eigen::VectorXd> active_modes;
  std::vector<double> active_omegas;
  active_modes.reserve(nmodes);
  active_omegas.reserve(nmodes);

  for (Index s = 0; s < nmodes; s++) {
    if (all_norms[s] > tol) {
      active_modes.push_back(std::move(all_modes[s]));
      active_omegas.push_back(omegas(s));
    }
  }

  screening_omegas_ = Eigen::VectorXd::Zero(Index(active_omegas.size()));
  for (Index s = 0; s < screening_omegas_.size(); s++) {
    screening_omegas_(s) = active_omegas[std::size_t(s)];
  }

  screening_modes_ = std::move(active_modes);
  screening_cached_ = true;
}

void RPA_UKS::ShiftUncorrectedEnergies(Eigen::VectorXd& energies,
                                       const Eigen::VectorXd& dftenergies,
                                       Index homo, Index qpmin, Index gwsize) {
  out_of_window_.Apply(energies, dftenergies, rpamin_, rpamax_, homo, qpmin,
                       qpmin + gwsize - 1);
}
Eigen::MatrixXd RPA_UKS::ResponseSum(const SpinWeightsFn& weights) const {
  // The dielectric matrix lives in the auxiliary basis; alpha and beta
  // channels add, each hole slice v contributing M_v^T diag(w_v) M_v.
  const Index size = Mmn_.alpha.auxsize();
  const Index n_occ_alpha = std::max<Index>(0, homo_alpha_ + 1 - rpamin_);
  const Index n_occ_beta = std::max<Index>(0, homo_beta_ + 1 - rpamin_);
  const Index n_unocc_alpha = std::max<Index>(0, rpamax_ - homo_alpha_);
  const Index n_unocc_beta = std::max<Index>(0, rpamax_ - homo_beta_);
  const Index nalpha = (n_unocc_alpha > 0) ? n_occ_alpha : 0;
  const Index nbeta = (n_unocc_beta > 0) ? n_occ_beta : 0;
  OptionalTiming timing(timings_, "screening: epsilon (RPA)");

  auto spin_of = [&](Index block, Index& m) {
    if (block < nalpha) {
      m = block;
      return false;
    }
    m = block - nalpha;
    return true;
  };

  return WeightedGram::Compute(
      nalpha + nbeta, size,
      [&](Index block, Eigen::VectorXd& w) {
        Index m = 0;
        const bool beta = spin_of(block, m);
        weights(beta, m, w);
      },
      [&](Index block, const WeightedGram::RowsVisitor& use) {
        Index m = 0;
        const bool beta = spin_of(block, m);
        if (beta) {
          use(Mmn_.beta[m].bottomRows(n_unocc_beta));
        } else {
          use(Mmn_.alpha[m].bottomRows(n_unocc_alpha));
        }
      });
}

template <bool imag>
Eigen::MatrixXd RPA_UKS::calculate_epsilon(double frequency) const {
  const double freq2 = frequency * frequency;
  const double eta2 = eta_ * eta_;
  Eigen::MatrixXd result =
      ResponseSum([&](bool beta, Index m_level, Eigen::VectorXd& denom) {
        const Eigen::VectorXd& energies =
            beta ? energies_beta_ : energies_alpha_;
        const Index homo = beta ? homo_beta_ : homo_alpha_;
        const Index n_unocc = rpamax_ - homo;
        // Bare particle-hole energy differences Delta_e = e_a - e_i
        const Eigen::ArrayXd deltaE =
            energies.tail(n_unocc).array() - energies(m_level);
        if (imag) {
          // 2 Delta_e / (Delta_e^2 + w^2) per spin channel (the closed-shell
          // code has 4: both spins folded into one channel)
          denom = 2.0 * deltaE / (deltaE.square() + freq2);
        } else {
          // (Delta_e - w)/((Delta_e - w)^2 + eta^2)
          //   + (Delta_e + w)/((Delta_e + w)^2 + eta^2)
          Eigen::ArrayXd deltEf = deltaE - frequency;
          Eigen::ArrayXd sum = deltEf / (deltEf.square() + eta2);
          deltEf = deltaE + frequency;
          sum += deltEf / (deltEf.square() + eta2);
          denom = sum;
        }
      });
  // epsilon = 1 - v chi0: the identity on the auxiliary-space diagonal
  result.diagonal().array() += 1.0;
  return result;
}

template Eigen::MatrixXd RPA_UKS::calculate_epsilon<true>(
    double frequency) const;
template Eigen::MatrixXd RPA_UKS::calculate_epsilon<false>(
    double frequency) const;

Eigen::MatrixXd RPA_UKS::calculate_epsilon_r(
    std::complex<double> frequency) const {
  // Real part at complex frequency z = w + i*gamma: terms of
  // Re[1 / (w + i*gamma - Delta_e)] and its partner at -Delta_e.
  const double sigma_1 = std::pow(frequency.imag() + eta_, 2);
  const double sigma_2 = std::pow(frequency.imag() - eta_, 2);
  Eigen::MatrixXd result =
      ResponseSum([&](bool beta, Index m_level, Eigen::VectorXd& chi) {
        const Eigen::VectorXd& energies =
            beta ? energies_beta_ : energies_alpha_;
        const Index homo = beta ? homo_beta_ : homo_alpha_;
        const Index n_unocc = rpamax_ - homo;
        const Eigen::ArrayXd deltaE =
            energies.tail(n_unocc).array() - energies(m_level);
        const Eigen::ArrayXd deltaEm = frequency.real() - deltaE;
        const Eigen::ArrayXd deltaEp = frequency.real() + deltaE;
        // the overall factor -1 (restricted convention per spin channel)
        // is folded into the weights
        chi = -(deltaEm * (deltaEm.abs2() + sigma_1).inverse() -
                deltaEp * (deltaEp.abs2() + sigma_2).inverse());
      });
  result.diagonal().array() += 1.0;
  return result;
}

const RPA_UKS::rpa_eigensolution& RPA_UKS::Diagonalize_H2p() const {
  if (h2p_cached_) {
    return h2p_solution_cache_;
  }
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

  const SymmetricEigenSystem es = Diagonalize_H2p_C(C);

  sol.omega = Eigen::VectorXd::Zero(es.values.size());
  sol.omega = es.values.cwiseSqrt();
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
        Omega_sqrt_inv(s) * AmB_sqrt.cwiseProduct(es.vectors.col(s));
  }

  h2p_solution_cache_ = std::move(sol);
  h2p_cached_ = true;
  return h2p_solution_cache_;
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
          kUKSApBPrefactor *
          Mmn_.alpha[v1].middleRows(n_occ_alpha, n_unocc_alpha) * Mmn_v2T;
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
          kUKSApBPrefactor *
          Mmn_.beta[v1].middleRows(n_occ_beta, n_unocc_beta) * Mmn_v2T;
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
          kUKSApBPrefactor *
          Mmn_.alpha[v_alpha].middleRows(n_occ_alpha, n_unocc_alpha) *
          Mmn_beta_T;
    }
  }

  // Symmetrize alpha-alpha block from its computed lower triangle.
  ApB.block(0, 0, size_alpha, size_alpha)
      .template triangularView<Eigen::StrictlyUpper>() =
      ApB.block(0, 0, size_alpha, size_alpha)
          .transpose()
          .template triangularView<Eigen::StrictlyUpper>();

  // Symmetrize beta-beta block from its computed lower triangle.
  ApB.block(size_alpha, size_alpha, size_beta, size_beta)
      .template triangularView<Eigen::StrictlyUpper>() =
      ApB.block(size_alpha, size_alpha, size_beta, size_beta)
          .transpose()
          .template triangularView<Eigen::StrictlyUpper>();

  // Copy the mixed alpha-beta block into the beta-alpha block.
  ApB.block(size_alpha, 0, size_beta, size_alpha) =
      ApB.block(0, size_alpha, size_alpha, size_beta).transpose();

  // Add the diagonal bare transition energies.
  ApB.diagonal() += Calculate_H2p_AmB();

  return ApB;
}

SymmetricEigenSystem RPA_UKS::Diagonalize_H2p_C(
    const Eigen::MatrixXd& C) const {
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Diagonalizing two-particle Hamiltonian "
      << std::flush;

  const SymmetricEigenSystem es = SymmetricEigen(C);

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Diagonalization done " << std::flush;

  const double minCoeff = es.values.minCoeff();
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