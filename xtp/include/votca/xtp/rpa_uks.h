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

 #pragma once
#ifndef VOTCA_XTP_RPA_UKS_H
#define VOTCA_XTP_RPA_UKS_H

#include <complex>
#include <vector>

#include "eigen.h"
#include "logger.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

/**
 * \brief Unrestricted RPA helper for spin-resolved GW screening.
 *
 * This class is the unrestricted/open-shell analogue of the existing restricted
 * RPA class. It constructs the independent-particle polarizability and related
 * dielectric quantities from explicit alpha and beta particle-hole channels.
 *
 * Conceptually, for a collinear unrestricted reference,
 *
 *   chi0(w) = chi0_alpha(w) + chi0_beta(w)
 *
 * where each spin contribution contains only same-spin particle-hole
 * excitations:
 *
 *   i_alpha -> a_alpha
 *   i_beta  -> a_beta
 *
 * The Coulomb interaction is spin-independent, so the resulting dielectric
 * matrix is built from the sum of both spin channels. This is exactly what is
 * needed later for unrestricted GW, where Sigma^alpha and Sigma^beta share the
 * same screened interaction W but differ in exchange and in the external
 * spin-resolved states on which Sigma acts.
 *
 * Compared to the restricted closed-shell implementation:
 *
 * - we do not assume one doubly occupied spatial-orbital ladder,
 * - we do not insert a blanket closed-shell spin-degeneracy factor,
 * - instead, alpha and beta channels are treated explicitly and summed.
 *
 * The class also provides a two-particle Hamiltonian construction
 * (Diagonalize_H2p), mainly mirroring the restricted implementation. For the
 * initial unrestricted GW work, the dielectric functions are the more critical
 * part; the H2p diagonalization is included for completeness and future use.
 */
class RPA_UKS {
 public:
  /**
   * \param log  Logger used for diagnostic output
   * \param Mmn  Spin-resolved three-center integrals in the MO basis
   *
   * Mmn_.alpha and Mmn_.beta contain the RI three-center objects transformed
   * with alpha and beta MOs, respectively.
   */
  RPA_UKS(Logger& log, const TCMatrix_gwbse_spin& Mmn) : log_(log), Mmn_(Mmn) {}

  /**
   * \brief Configure orbital window and spin-resolved HOMO indices.
   *
   * The same global orbital window [rpamin_, rpamax_] is used for both spin
   * channels, but the occupied/virtual split is spin dependent:
   *
   *   alpha: occupied <= homo_alpha_, virtual >= homo_alpha_ + 1
   *   beta : occupied <= homo_beta_,  virtual >= homo_beta_  + 1
   *
   * This keeps the interface close to the restricted implementation while still
   * allowing different alpha and beta occupations.
   */
  void configure(Index homo_alpha, Index homo_beta, Index rpamin,
                 Index rpamax) {
    homo_alpha_ = homo_alpha;
    homo_beta_ = homo_beta;
    rpamin_ = rpamin;
    rpamax_ = rpamax;
  }

  /**
   * \brief Small positive broadening used in the real-frequency response.
   */
  double getEta() const { return eta_; }

  /**
   * \brief Dielectric matrix on the imaginary frequency axis.
   *
   * This evaluates the analogue of the restricted epsilon(iw), but with
   * explicit alpha and beta particle-hole sums.
   */
  Eigen::MatrixXd calculate_epsilon_i(double frequency) const {
    return calculate_epsilon<true>(frequency);
  }

  /**
   * \brief Dielectric matrix on the real frequency axis.
   *
   * The implementation uses the same algebraic structure as the restricted
   * version, but accumulates alpha and beta channel contributions explicitly.
   */
  Eigen::MatrixXd calculate_epsilon_r(double frequency) const {
    return calculate_epsilon<false>(frequency);
  }

  /**
   * \brief Real part of the dielectric matrix at complex frequency.
   *
   * This is used in the same contexts as in the restricted code path, but with
   * unrestricted spin resolution.
   */
  Eigen::MatrixXd calculate_epsilon_r(std::complex<double> frequency) const;

  /**
   * \brief Access the current alpha RPA input energies used in the response.
   *
   * These are typically DFT orbital energies corrected by available GW energies
   * in the qp window, with outer states shifted by a constant maximal
   * correction, matching the philosophy of the restricted implementation.
   */
  const Eigen::VectorXd& getRPAInputEnergiesAlpha() const {
    return energies_alpha_;
  }

  /**
   * \brief Access the current beta RPA input energies used in the response.
   */
  const Eigen::VectorXd& getRPAInputEnergiesBeta() const {
    return energies_beta_;
  }

  /**
   * \brief Set alpha and beta RPA input energies directly.
   *
   * This is the unrestricted counterpart of setting one restricted RPA energy
   * ladder. The vectors are assumed to span the configured range
   * [rpamin_, rpamax_].
   */
  void setRPAInputEnergies(const Eigen::VectorXd& rpaenergies_alpha,
                           const Eigen::VectorXd& rpaenergies_beta) {
    energies_alpha_ = rpaenergies_alpha;
    energies_beta_ = rpaenergies_beta;
  }

  /**
   * \brief Update alpha and beta RPA input energies from DFT and GW energies.
   *
   * Inside the qp window, available GW quasiparticle energies replace the DFT
   * ones. Outside that window, states are shifted by a constant correction
   * equal to the largest correction seen in the occupied or virtual GW region,
   * respectively. This mirrors the existing restricted strategy.
   */
  void UpdateRPAInputEnergies(const Eigen::VectorXd& dftenergies_alpha,
                              const Eigen::VectorXd& dftenergies_beta,
                              const Eigen::VectorXd& gwaenergies_alpha,
                              const Eigen::VectorXd& gwaenergies_beta,
                              Index qpmin);

  /**
   * \brief Solution of the two-particle Hamiltonian in the unrestricted basis.
   *
   * The particle-hole basis is ordered as
   *
   *   [ all alpha excitations | all beta excitations ]
   *
   * with alpha excitations enumerated first and beta excitations appended
   * afterwards. This ordering is used consistently in Calculate_H2p_AmB() and
   * Calculate_H2p_ApB().
   */
  struct rpa_eigensolution {
    Eigen::VectorXd omega;
    Eigen::MatrixXd XpY;
    double ERPA_correlation;
  };

  /**
   * \brief Diagonalize the unrestricted two-particle Hamiltonian.
   *
   * This constructs the unrestricted analogue of the restricted H2p problem.
   * The diagonal AmB part contains bare particle-hole energy differences, while
   * ApB adds Coulomb coupling between particle-hole excitations from both spin
   * channels.
   */
  rpa_eigensolution Diagonalize_H2p() const;

 private:
  /**
   * \brief Internal implementation of epsilon(w).
   *
   * imag=true  -> imaginary frequency branch
   * imag=false -> real frequency branch
   *
   * The result is accumulated as a sum of alpha and beta independent-particle
   * contributions in the auxiliary basis.
   */
  template <bool imag>
  Eigen::MatrixXd calculate_epsilon(double frequency) const;

  /**
   * \brief Construct the diagonal AmB block of the two-particle Hamiltonian.
   *
   * For each particle-hole excitation v -> c, this contributes the bare
   * excitation energy
   *
   *   Delta_e = e_c - e_v
   *
   * evaluated separately for alpha and beta channels and then concatenated into
   * one vector.
   */
  Eigen::VectorXd Calculate_H2p_AmB() const;

  /**
   * \brief Construct the ApB block of the unrestricted two-particle Hamiltonian.
   *
   * This contains Coulomb coupling between particle-hole excitations. In the
   * unrestricted case, all four spin blocks are present:
   *
   *   alpha-alpha, alpha-beta, beta-alpha, beta-beta
   *
   * because the Coulomb interaction couples total density fluctuations and is
   * spin independent.
   */
  Eigen::MatrixXd Calculate_H2p_ApB() const;

  /**
   * \brief Diagonalize the symmetrized two-particle matrix C.
   *
   * This helper mirrors the restricted implementation and performs the final
   * Hermitian diagonalization after symmetrization with sqrt(AmB).
   */
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> Diagonalize_H2p_C(
      const Eigen::MatrixXd& C) const;

  /**
   * \brief Shift uncorrected states outside the GW window.
   *
   * Occupied states below qpmin are shifted by the largest occupied correction
   * found within the GW window, and virtual states above qpmax are shifted by
   * the largest virtual correction. This is done independently for each spin
   * channel.
   */
  void ShiftUncorrectedEnergies(Eigen::VectorXd& energies,
                                const Eigen::VectorXd& dftenergies, Index homo,
                                Index qpmin, Index gwsize);

  /**
   * \brief Return the largest absolute GW correction in a given index range.
   */
  double getMaxCorrection(const Eigen::VectorXd& energies,
                          const Eigen::VectorXd& dftenergies, Index min,
                          Index max) const;

  // Spin-resolved HOMO indices in the full orbital numbering
  Index homo_alpha_;
  Index homo_beta_;

  // Shared orbital window used to build the RPA quantities
  Index rpamin_;
  Index rpamax_;

  // Small broadening parameter used in real-frequency expressions
  const double eta_ = 0.0001;

  // Current alpha/beta energies used in the response over [rpamin_, rpamax_]
  Eigen::VectorXd energies_alpha_;
  Eigen::VectorXd energies_beta_;

  Logger& log_;

  // Spin-resolved three-center integrals:
  // Mmn_.alpha for alpha MOs, Mmn_.beta for beta MOs
  const TCMatrix_gwbse_spin& Mmn_;
};

}  // namespace xtp
}  // namespace votca

#endif