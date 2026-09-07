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
#ifndef VOTCA_XTP_EWALDREALSPACESUM_H
#define VOTCA_XTP_EWALDREALSPACESUM_H

// Standard includes
#include <unordered_map>
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "ewaldrealspaceinteractor.h"
#include "ewaldregistry.h"

/**
 * \brief Real-space term of a 3D Ewald sum: the total erfc(alpha*r)-screened
 *        field at a target position from every registered segment, summed
 *        over periodic images, converged by adaptive radial shells.
 *
 * This is a fresh implementation of the algorithm identified in the legacy
 * xtp/ewald/ code's RThread::FP_FieldCalc (real-space contribution grows
 * shell-by-shell until the per-shell contribution becomes negligible), not
 * a port of it. In particular the convergence criterion is re-derived: the
 * legacy code's own "test dipole" / RMS-count construction was never fully
 * understood even by the people who wrote it, so this class instead tracks
 * the plain magnitude of each shell's own field contribution and stops once
 * that drops below an absolute tolerance (with a minimum radius enforced
 * throughout, so an accidentally-small early shell can't trigger early
 * stopping). This is a considered simplification, not a proven-equivalent
 * reformulation -- see the project notes on validating its own convergence
 * behaviour (does the total field actually stabilize as the tolerance
 * tightens?) before trusting it on a production system.
 *
 * Scope: this class computes the *intermolecular* real-space field only.
 * Interactions between two sites in the same registered segment (whether
 * at zero translation, i.e. genuinely intramolecular, or at a nonzero
 * translation, i.e. that segment's own periodic image) are excluded
 * entirely and left to whatever already handles intramolecular coupling
 * -- this mirrors the existing split in eeInteractor between
 * ApplyStaticField/ApplyInducedField (intermolecular) and
 * Cholesky_IntraSegment/FillTholeInteraction (intramolecular), rather than
 * inventing a new boundary.
 *
 * Performance: AddFieldAt caches the real, geometry-determined neighbor
 * set (which (source segment, periodic-translation) pairs actually
 * contributed to a given target's own converged sum) per (target site,
 * source charge state) pair, keyed by that site's own address together
 * with source_state. The first call for a given (target, state) pair
 * runs the full shell-by-shell convergence search (scanning every
 * registered segment each shell); every subsequent call for that SAME
 * pair reuses the cached (source, translation) list directly, skipping
 * the search entirely. This mirrors legacy PolarBackground's own RThread::
 * PolarNbs() mechanism (built once, at its own first SOR iteration, then
 * reused for every later one) rather than being a novel optimization --
 * a direct comparison against a legacy log for the same real system
 * found this class's own repeated full-system rescan, once per PCG
 * iteration, to be the dominant cause of a real ~8x wall-clock slowdown
 * relative to legacy's own SOR solve at a comparable iteration count
 * (with thread count and iteration count both directly ruled out as
 * explanations first). This caching is safe specifically because site
 * POSITIONS are fixed for the lifetime of a single PCG solve (only
 * induced dipole VALUES change between iterations, which this cache
 * never stores) -- it would NOT be safe to reuse across a call sequence
 * where target or source positions themselves change.
 *
 * Units: bohr, Hartree atomic units throughout, matching PolarSite and
 * EwaldRealSpaceInteractor.
 */

namespace votca {
namespace xtp {

class EwaldRealSpaceSum {
 public:
  // box: columns are the lattice vectors a, b, c (matching
  //   Topology::getBox()'s own convention).
  // registry: source of every segment's multipole/induced-dipole state.
  //   Stored by reference; the registry must outlive this object, and
  //   FieldAt() reads whatever is currently in it (e.g. updated induced
  //   dipoles from a prior SCF iteration).
  // alpha, thole_a: passed through to the underlying
  //   EwaldRealSpaceInteractor.
  // r_min: minimum shell radius (bohr) before convergence may be declared,
  //   regardless of how small any individual shell's contribution is.
  // field_tol: convergence threshold (atomic units) on a shell's own
  //   contribution to the field magnitude.
  // shell_width: radial bin width (bohr) used to group periodic images
  //   into shells.
  // n_max: safety cap on the lattice-translation search box (translations
  //   with |na|,|nb|,|nc| <= n_max are considered candidates). If
  //   convergence is not reached within this box, FieldAt() throws rather
  //   than silently returning an unconverged result.
  EwaldRealSpaceSum(const Eigen::Matrix3d& box, const EwaldRegistry& registry,
                    double alpha, double thole_a, double r_min,
                    double field_tol, double shell_width = 0.945,
                    Index n_max = 15);

  // Accumulates the total intermolecular real-space field into target's
  // own V()/V_noE() accumulators (via EwaldRealSpaceInteractor, matching
  // its convention), generated by every segment in the registry at charge
  // state source_state, excluding target_segment_id itself (see class
  // documentation on scope). target's position (target.getPos()) is used
  // as the field evaluation point; target itself need not already be
  // registered in registry_. Throws std::runtime_error if convergence is
  // not reached within n_max shells.
  template <enum Estatic CE>
  void AddFieldAt(Index target_segment_id, PolarSite& target,
                  EwaldChargeState source_state) const;

 private:
  // One periodic image translation vector, tagged with its distance from
  // the origin so shells can be built by sorting once.
  struct Translation {
    Eigen::Vector3d t;
    double r;
  };
  std::vector<Translation> GenerateSortedTranslations() const;

  // std::pair has no default std::hash specialization; EwaldChargeState
  // (an enum class) does, via std::underlying_type, since C++14, so this
  // only needs to combine the two.
  struct PairHash {
    std::size_t operator()(
        const std::pair<const PolarSite*, EwaldChargeState>& key) const {
      return std::hash<const PolarSite*>()(key.first) ^
            (std::hash<EwaldChargeState>()(key.second) << 1);
    }
  };

  Eigen::Matrix3d box_;
  const EwaldRegistry& registry_;
  EwaldRealSpaceInteractor interactor_;
  double r_min_;
  double field_tol_;
  double shell_width_;
  Index n_max_;

  // Generated once at construction (depends only on the box, not on the
  // target position), reused by every FieldAt() call.
  std::vector<Translation> translations_;

  // Every (source_id, translation_idx) pair that genuinely contributed to
  // a given target's own converged sum -- i.e. the real, geometry-
  // determined neighbor list, analogous to legacy PolarBackground's own
  // RThread::PolarNbs(). Built once per distinct (target site, source
  // charge state) pair (keyed by the site's own address plus
  // source_state -- the address alone is NOT enough: this class's
  // registered segments could in principle be queried at a different
  // EwaldChargeState for the same physical target site, and geometric
  // neighbors depend on which segments are actually registered at that
  // state via Has(), not on the target's own identity alone), then
  // reused on every subsequent call for that same (target, state) pair.
  // The motivating caller (EwaldPeriodicDipoleOperator, via its own
  // targets_ vector built once and reused every RawMultiply call) always
  // uses EwaldChargeState::Neutral, so this distinction is not
  // exercised by anything in this codebase today, but the key is chosen
  // to be genuinely correct rather than correct only for the one
  // access pattern that happens to exist right now.
  mutable std::unordered_map<std::pair<const PolarSite*, EwaldChargeState>,
                             std::vector<std::pair<Index, Index>>,
                             PairHash>
      neighbor_cache_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDREALSPACESUM_H
