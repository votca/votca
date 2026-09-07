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
#ifndef VOTCA_XTP_EWALDRECIPROCALSPACESUM_H
#define VOTCA_XTP_EWALDRECIPROCALSPACESUM_H

// Standard includes
#include <complex>
#include <functional>
#include <vector>

// Local VOTCA includes
#include "eeinteractor.h"
#include "eigen.h"
#include "ewaldregistry.h"

/**
 * \brief Reciprocal-space term of a 3D Ewald sum: the field at a target
 *        position from the Gaussian-screened periodic image of every
 *        registered charge and (static or induced) dipole.
 *
 * This is derived from scratch, not ported from the legacy xtp/ewald/ code.
 * Ewald's trick surrounds each point multipole with a canceling Gaussian
 * charge distribution; that periodic Gaussian density's own potential
 * satisfies Poisson's equation in Fourier space (phi_k = 4*pi*rho_k/k^2),
 * which gives, for a charge+dipole source set:
 *
 *   S(k)   = sum_j [ q_j - i*(k.mu_j) ] * exp(-i*k.r_j)      (structure
 *                                                              factor)
 *   phi(r) = (4*pi/V) * sum_{k!=0} (1/k^2) * exp(-k^2/4a^2) * S(k) *
 *            exp(i*k.r)
 *   E(r)   = (4*pi/V) * sum_{k!=0} (k/k^2) * exp(-k^2/4a^2) *
 *            Im[ S(k) * exp(i*k.r) ]
 *
 * summed over the *full* set of k!=0 (not a k/-k half-space with an extra
 * factor of 2) -- pairing k and -k terms is exactly what makes this sum
 * real by construction, and is worth stating explicitly since getting it
 * wrong (e.g. by introducing a spurious factor of 2, or its inverse) would
 * not be caught by compilation, only by a numerical check against a known
 * reference. See the project notes on validating this class for exactly
 * that reason.
 *
 * k-vector truncation is a simple |k| < k_max spherical cutoff, not the
 * legacy code's per-vector "graded" adaptive truncation. That grading is
 * self-contained (depends only on structure-factor magnitude, not on any
 * real-space convergence state -- confirmed by tracing the legacy code
 * directly), so it is addable later as a genuinely separable enhancement
 * to the truncation strategy, on top of this simple-cutoff implementation,
 * without needing to revisit the underlying k-vector definitions or field
 * formula above.
 *
 * Scope: unlike EwaldRealSpaceSum, this includes EVERY registered site --
 * no exclusion of the target's own segment at all, even for the
 * intramolecular case. This was not always true of this class: an
 * earlier version excluded the target's own segment entirely, modeled on
 * EwaldRealSpaceSum's own real-space exclusion, based on a mistaken
 * assumption that the two sums should treat exclusion the same way.
 * Tracing legacy PolarBackground's own reciprocal-space code directly
 * (PolarBackground::KThread::SP_SFactorCalc/FP_KFieldCalc for the
 * permanent case, SU_SFactorCalc/FU_KFieldCalc for the induced case)
 * shows it never excludes anything from its own global structure factor
 * -- every registered site contributes, and the resulting field is
 * applied back to every site with no segment-awareness at all. This
 * matters specifically because of how the Ewald split works: erfc(a*r)/r
 * + [this class's own reciprocal-space contribution] = 1/r EXACTLY (no
 * approximation) -- so if EwaldRealSpaceSum's own real-space sum
 * legitimately excludes a pair (as it does for same-segment,
 * zero-translation pairs), the reciprocal-space contribution to that
 * SAME pair must still be included here for the two Ewald halves to
 * combine correctly wherever they're meant to (see
 * EwaldPeriodicDipoleOperator's own AddIntraSegmentCoupling for the
 * intramolecular case this enables, and its own class documentation for
 * the fuller account of how this was discovered).
 *
 * Performance note: S(k) is a sum over every registered site, and is
 * identical for every target sharing the same source_state -- no
 * per-target correction of any kind is needed now (there is nothing
 * target-specific left to correct for), so AddFieldAtMany's own benefit
 * over repeated AddFieldAt calls is now simply computing S(k) once
 * per k-vector and reusing it across every target in the batch, rather
 * than once per (k-vector, target) pair. This is an O(N_k * N) batch,
 * rather than the O(N_targets * N_k * N) that repeated single-target
 * AddFieldAt calls would cost -- the difference that matters once this
 * is driven from an iterative solve calling it once per site per
 * iteration. AddFieldAt (single target) is kept for convenience/testing
 * and simply delegates to AddFieldAtMany with a one-element batch; it
 * does not itself get the batching benefit, since there's nothing to
 * share it with.
 *
 * Units: bohr, Hartree atomic units throughout.
 */

namespace votca {
namespace xtp {

class EwaldReciprocalSpaceSum {
 public:
  // box: columns are the lattice vectors a, b, c.
  // registry: source of every segment's multipole/induced-dipole state.
  //   Stored by reference; must outlive this object.
  // alpha: Ewald splitting parameter, shared with whatever
  //   EwaldRealSpaceInteractor/EwaldRealSpaceSum instance is used
  //   alongside this class for the same background.
  // k_max: spherical cutoff (bohr^-1) on |k|.
  EwaldReciprocalSpaceSum(const Eigen::Matrix3d& box,
                          const EwaldRegistry& registry, double alpha,
                          double k_max);

  // Number of k-vectors this instance will sum over (fixed at
  // construction) -- useful for logging/progress reporting before or
  // during a call to AddFieldAt/AddFieldAtMany, since this count (not
  // k_max itself) is what actually determines the cost of a call: it
  // scales with the cube of k_max relative to the box's own reciprocal
  // lattice spacing, so a k_max that looks modest can still generate an
  // unexpectedly large count for a large box.
  std::size_t NumKVectors() const { return kvectors_.size(); }

  // The position-independent 3x3 matrix M such that a site's own dipole
  // moment mu produces a spurious reciprocal-space "self-field" -M*mu at
  // its own position -- a known Ewald-summation artifact: since this
  // class's own structure factor S(k) never excludes any site, not even
  // the one a field is being evaluated at (see class documentation), a
  // site with a nonzero dipole moment feels a nonzero contribution from
  // itself. This does NOT happen for a pure charge (monopole) term --
  // that self-contribution is exactly zero by construction (a charge's
  // own term in S(k), evaluated back at its own position, is purely real
  // -- q*exp(-i*k.r)*exp(i*k.r) = q -- and this class's own field formula
  // only ever keeps the imaginary part) -- so this matrix, and the
  // correction it enables, only ever matters for dipoles (static or
  // induced), never for bare charges.
  //
  // M_ab = (4*pi/V) * sum_{k!=0} (k_a*k_b/k^2) * exp(-k^2/4*alpha^2)
  //
  // -- the same k-vector set, weight, and prefactor AddFieldAtMany itself
  // uses, so M is exactly the numerical self-term this class's own sum
  // actually produces (not an analytic/closed-form approximation to it,
  // which would generally differ once k_max is finite rather than
  // infinite). Symmetric by construction (k_a*k_b = k_b*k_a). Neither
  // this class nor EwaldRealSpaceSum currently subtracts this
  // contribution automatically anywhere -- confirmed (by tracing its own
  // reciprocal-space code directly) that legacy PolarBackground does not
  // either, so the two remain consistent with each other as they
  // currently stand. A caller wanting the physically-corrected (self-
  // interaction-free) field must subtract SelfFieldMatrix()*mu itself,
  // for whichever site's own dipole moment mu is being evaluated.
  Eigen::Matrix3d SelfFieldMatrix() const;

  // Progress callback signature: called as (k_vectors_done,
  // k_vectors_total) periodically (not necessarily every single
  // k-vector) during the expensive part of AddFieldAtMany. Pass an empty
  // std::function (the default) for no callback at all.
  using ProgressCallback = std::function<void(std::size_t, std::size_t)>;

  // Accumulates the total reciprocal-space field into target's own
  // V()/V_noE() accumulators (matching EwaldRealSpaceSum's convention),
  // generated by every registered segment at charge state source_state
  // (see class documentation on scope -- no segment is excluded, not
  // even target's own). Convenience wrapper around AddFieldAtMany for a
  // single target; prefer AddFieldAtMany when computing the field at
  // many sites at once (see the class documentation's performance note).
  template <enum Estatic CE>
  void AddFieldAt(PolarSite& target, EwaldChargeState source_state) const;

  // Batched form of AddFieldAt: computes the field at every target site
  // in targets, sharing a single pass over every k-vector's structure
  // factor across the whole batch. target pointers must be non-null and
  // must remain valid for the duration of the call; the same site may
  // appear more than once (field is simply accumulated twice in that
  // case, same as two separate AddFieldAt calls would do). progress, if
  // non-empty, is called periodically during the k-vector loop (see
  // ProgressCallback).
  template <enum Estatic CE>
  void AddFieldAtMany(const std::vector<PolarSite*>& targets,
                      EwaldChargeState source_state,
                      const ProgressCallback& progress = ProgressCallback()) const;

 private:
  struct KVector {
    Eigen::Vector3d k;
    double k2;
  };
  std::vector<KVector> GenerateKVectors() const;

  // Structure factor S(k) (see class documentation), one entry per
  // kvectors_, summed over every site of every registered segment at
  // source_state. This is the expensive O(N_k * N_sites) loop
  // AddFieldAtMany's own progress callback reports on.
  std::vector<std::complex<double>> TotalStructureFactors(
      EwaldChargeState source_state, const ProgressCallback& progress) const;

  Eigen::Matrix3d box_;
  double volume_;
  const EwaldRegistry& registry_;
  double alpha_;
  double k_max_;

  std::vector<KVector> kvectors_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDRECIPROCALSPACESUM_H
