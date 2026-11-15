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

#pragma once
#ifndef VOTCA_XTP_INACTIVESCREENING_H
#define VOTCA_XTP_INACTIVESCREENING_H

// Standard includes
#include <string>
#include <vector>

// Local VOTCA includes
#include "aobasis.h"
#include "eigen.h"
#include "orbitals.h"
#include "threecenter.h"

namespace votca {
namespace xtp {

/**
 * \brief Occupied orbitals of the inactive region of a projection-based
 * embedding, with the virtual space they are excited into.
 *
 * The occupied orbitals are the inactive part of the occupied space of the
 * full reference calculation, brought into semicanonical form (the Fock
 * operator is diagonal within them). The embedded active orbitals are the
 * semicanonical orbitals of the complement (eigenfunctions of Q F_full Q),
 * so F_full still couples active and inactive occupied orbitals, and the
 * response summed over the two blocks separately differs from the full one
 * at first order in that coupling. The virtual orbitals are the full
 * reference virtuals, equal to the embedded virtuals at convergence.
 */
struct InactiveSpace {
  Eigen::MatrixXd occupied;  // AO x n_inactive
  Eigen::VectorXd occupied_energies;
  Eigen::MatrixXd virtuals;  // AO x n_virtual
  Eigen::VectorXd virtual_energies;
  // max |Q^2 - Q| of the projector of the inactive density onto the
  // occupied space; 0 if the inactive density is made of occupied orbitals
  double idempotency_error = 0.0;

  Index size() const { return occupied.cols(); }
  bool empty() const { return occupied.cols() == 0; }
};

/**
 * \brief What GWBSE needs to add the static screening of the inactive
 * region: the mode and the orbitals of the full reference calculation that
 * enter the response (see InactiveScreening::Prepare).
 *
 *  - block: the semicanonical inactive occupied orbitals (InactiveSpace);
 *    Y = X_B. Misses the terms of first order in the Fock coupling between
 *    active and inactive occupied orbitals; cost proportional to the number
 *    of inactive orbitals.
 *  - exact: all canonical occupied orbitals of the full calculation;
 *    Y = X_full - X_A, with X_A from the embedded active orbitals. Gives the
 *    static screened interaction of the full calculation exactly; cost
 *    proportional to all occupied orbitals.
 */
struct InactiveScreeningInput {
  enum class Mode { None, Block, Exact };
  Mode mode = Mode::None;
  Eigen::MatrixXd occupied;
  Eigen::VectorXd occupied_energies;
  Eigen::MatrixXd virtuals;
  Eigen::VectorXd virtual_energies;
  Index n_inactive = 0;
  double idempotency_error = 0.0;
  // block mode with a spatial cutoff: inactive orbitals left out
  Index n_far = 0;
  double radius = 0.0;  // bohr; 0: all
  // complete G: the inactive occupied orbitals (semicanonical) also enter
  // the self-energy; they are put in front of the active ones
  // (InactiveScreening::CompleteGOrbitals) and excluded from the RPA
  bool complete_g = false;
  Eigen::MatrixXd inactive_occupied;
  Eigen::VectorXd inactive_energies;

  bool empty() const { return mode == Mode::None; }
  static Mode Parse(const std::string& name);
  static std::string Name(Mode mode);
};

/**
 * \brief Static screening by the inactive region of a PbE calculation.
 *
 * In the auxiliary frame of the stored three-centre integrals
 * (M = (mn|Q) T, T T^T = V^-1) the static RPA dielectric matrix is
 * 1 + X with X = sum_ia 4/(e_a - e_i) M_ia M_ia^T. Dressing the active
 * integrals with S = (1 + R)^(1/2), R = (1 + Y)^-1 - 1, for a positive
 * semidefinite Y gives
 *   S (1 + S X_A S)^-1 S = (1 + X_A + Y)^-1.
 * With Y = X_full - X_A ("exact") this is the static screened interaction
 * of the full calculation; with Y = X_B over the semicanonical inactive
 * orbitals ("block") it is that up to the active-inactive Fock coupling.
 * R enters GW and BSE like the reaction field of a polarizable environment
 * (GW::setReactionField, BSE::setReactionField), and combines with one:
 * 1 + R = ((1 + R_env)^-1 + Y)^-1, the Dyson equation in which the
 * inactive region and the environment polarize each other.
 */
class InactiveScreening {
 public:
  /**
   * \brief The inactive space from the full reference orbitals and the
   * inactive density stored by DFTEngine::EvaluateActiveRegion.
   *
   * `full` must still hold the MOs of the full reference calculation (not
   * the embedded ones). The lowest `ncore` semicanonical inactive orbitals
   * are left out (frozen core).
   */
  static InactiveSpace BuildSpace(const Orbitals& full, Index ncore = 0);

  /**
   * \brief X = sum_ia 4/(e_a - e_i) M_ia M_ia^T over the given occupied
   * and virtual orbitals, in the frame TCMatrix_gwbse::Fill stores its
   * integrals in for `auxbasis`.
   *
   * The occupied orbitals are processed in chunks of at most
   * `chunk_bytes` of integrals, so the full tensor is never stored.
   * With `imaginary_frequency` w > 0 the weights are 4D/(D^2 + w^2), the
   * response at i w as RPA::calculate_epsilon_i builds it.
   */
  static Eigen::MatrixXd Response(
      const AOBasis& auxbasis, const AOBasis& dftbasis,
      const Eigen::MatrixXd& occupied, const Eigen::VectorXd& occupied_energies,
      const Eigen::MatrixXd& virtuals, const Eigen::VectorXd& virtual_energies,
      double chunk_bytes = 4e9, double imaginary_frequency = 0.0);

  /// Response at several imaginary frequencies from one integral pass.
  /// `inv_sqrt`: the metric V^-1/2 of `auxbasis` (TCMatrix_gwbse::InvSqrt),
  /// if already known; otherwise computed once for all chunks.
  static std::vector<Eigen::MatrixXd> Responses(
      const AOBasis& auxbasis, const AOBasis& dftbasis,
      const Eigen::MatrixXd& occupied, const Eigen::VectorXd& occupied_energies,
      const Eigen::MatrixXd& virtuals, const Eigen::VectorXd& virtual_energies,
      double chunk_bytes, const std::vector<double>& imaginary_frequencies,
      const Eigen::MatrixXd& inv_sqrt = Eigen::MatrixXd());

  static Eigen::MatrixXd Response(const AOBasis& auxbasis,
                                  const AOBasis& dftbasis,
                                  const InactiveSpace& space,
                                  double chunk_bytes = 4e9,
                                  double imaginary_frequency = 0.0) {
    return Response(auxbasis, dftbasis, space.occupied, space.occupied_energies,
                    space.virtuals, space.virtual_energies, chunk_bytes,
                    imaginary_frequency);
  }

  /**
   * \brief The input for GWBSE, from the full reference orbitals (before
   * they are replaced by the embedded ones). Mode None returns an empty
   * input. Throws for a truncated basis.
   */
  static InactiveScreeningInput Prepare(const Orbitals& full,
                                        InactiveScreeningInput::Mode mode,
                                        bool complete_g = false,
                                        double radius = 0.0);

  /**
   * \brief Block mode with a spatial cutoff: the inactive LMOs whose charge
   * centroid lies within `radius` (bohr) of the centroid of an active LMO,
   * semicanonical among themselves (the Fock operator of the full
   * calculation diagonalized in their span). The response of the farther
   * inactive orbitals is left out; their coupling to the near ones enters
   * at the same order as the active-inactive coupling of block mode.
   * Needs the LMOs of the localization step.
   */
  static InactiveSpace NearSpace(const Orbitals& full, double radius,
                                 Index& n_far);

  /**
   * \brief The MOs GW-BSE works with for complete G: the inactive occupied
   * orbitals (semicanonical, KS energies), then the active occupied and the
   * virtual embedded MOs; the inactive levels the projector shifted to the
   * top are left out. All of them are occupied for the self-energy and
   * Vxc (full density); the RPA leaves out the first n_inactive holes.
   */
  static tools::EigenSystem CompleteGOrbitals(
      const tools::EigenSystem& embedded, Index n_active_occupied,
      const InactiveScreeningInput& input);

  /**
   * \brief Y, the response of the inactive region: X_B (block) or
   * X_full - X_A (exact), static or at an imaginary frequency. `embedded` holds
   * the embedded MOs as GW uses them (active occupied levels 0..homo); the
   * virtuals of X_A are the input's number of virtuals above the HOMO, which
   * leaves out the inactive levels the projector has shifted to the top.
   */
  static Eigen::MatrixXd Kernel(const AOBasis& auxbasis,
                                const AOBasis& dftbasis,
                                const InactiveScreeningInput& input,
                                const tools::EigenSystem& embedded, Index homo,
                                double chunk_bytes,
                                double imaginary_frequency = 0.0,
                                Index first_active = 0);

  /**
   * \brief Y at several imaginary frequencies, as GWBSE uses it: one
   * integral pass for the full (exact) or inactive (block) occupied
   * orbitals, with the full calculation's virtuals inside the RPA range
   * (homo+1..rpamax; the inactive levels the projector shifted above them
   * are left out). In exact mode X_A comes from `active`, the integrals GW
   * holds (rpa_energies for levels rpamin..rpamax, excluded holes left out
   * as in the RPA), so the active integrals are not built a second time. The
   * `virtual_cutoff` > 0 (Ha above the LUMO) keeps only the virtuals up to
   * it, in both terms of exact mode. The
   * static point is the imaginary frequency eta, as RPA::calculate_epsilon_r
   * at 0 (identical weights).
   */
  static std::vector<Eigen::MatrixXd> Kernels(
      const AOBasis& auxbasis, const AOBasis& dftbasis,
      const InactiveScreeningInput& input, const TCMatrix_gwbse& active,
      const Eigen::VectorXd& rpa_energies, Index homo, Index rpamin,
      Index rpamax, double chunk_bytes,
      const std::vector<double>& imaginary_frequencies,
      double virtual_cutoff = 0.0);

  /**
   * \brief Poles in the form of transitions that reproduce a response
   * at zero and at one imaginary frequency:
   *   Y(i w) = sum_k 4 W_k/(W_k^2 + w^2) c_k c_k^T,
   * exactly at w = 0 and w = w_i on the range of Y(0). With Y(0) = G G^T
   * (eigenvectors above `tolerance` times the largest eigenvalue) and
   * G^+ Y(i w_i) G^+T = Q diag(l) Q^T, l_k = W_k^2/(W_k^2 + w_i^2) and
   * c_k = sqrt(W_k / 4) (G Q)_k. For a genuine sum of transitions every l
   * lies in (0, 1); for a difference of two responses (exact mode) it may
   * not, and such modes keep only their static response (a pole at 1e3 Ha;
   * counted).
   */
  struct PoleFit {
    Eigen::VectorXd energies;
    Eigen::MatrixXd vectors;  // aux x poles
    Index dropped = 0;        // directions of Y(0) below the tolerance
    Index clipped = 0;        // modes with l outside (0, 1)
  };
  static PoleFit FitPoles(const Eigen::MatrixXd& Y0, const Eigen::MatrixXd& Yi,
                          double w_i, double tolerance = 1e-10);

  /**
   * \brief The reaction field R with 1 + R = ((1 + R_env)^-1 + Y)^-1;
   * without an environment (R_env empty) R = (1 + Y)^-1 - 1.
   */
  static Eigen::MatrixXd ReactionField(
      const Eigen::MatrixXd& Y,
      const Eigen::MatrixXd& R_env = Eigen::MatrixXd());
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_INACTIVESCREENING_H
