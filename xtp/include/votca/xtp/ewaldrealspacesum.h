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
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
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
  // screening_factor: sets the real-space distance cutoff to
  //   screening_factor / alpha. A pair separated by more than that
  //   contributes at most ~erfc(screening_factor) of the local field, so
  //   the default 6.0 truncates at ~2e-17 per pair -- far below any
  //   realistic field_tol, while discarding the large majority of
  //   (source segment, translation) pairs a periodic sum would otherwise
  //   evaluate (measured: 78.7% on a 1000-segment box, a 4.5x saving in
  //   real-space cost with the converged dipoles unchanged in every
  //   printed digit). Raise it to tighten the truncation at
  //   proportionally greater cost -- the number of surviving pairs grows
  //   as the cube -- or lower it to trade accuracy for speed
  //   deliberately. Values below ~4 start to matter at the 1e-11 level
  //   and should be checked against a reference before use.
  //
  // screening_factor is deliberately the LAST parameter rather than
  // sitting next to the other accuracy controls it belongs with: adding
  // it mid-signature silently reinterpreted the positional shell_width
  // and n_max arguments of existing callers, which still compiled and
  // then culled every pair. Keep new parameters at the end.
  //
  // foreground: (segment id, segment centroid position) pairs that are
  //   handled EXPLICITLY by a polar region rather than by this periodic
  //   sum, and so must be removed from it -- the "carving out" an MM/MM
  //   or QM/MM job performs. Each entry suppresses exactly the one
  //   periodic copy of that segment sitting at the given position; every
  //   OTHER lattice image of the same segment remains part of the
  //   background and is still summed.
  //
  //   Matched by position rather than by "is this the minimum image",
  //   deliberately. Two foreground segments can be separated by up to
  //   twice the region cutoff, so for a target near the edge of the
  //   foreground the minimum image of another foreground segment can be
  //   a DIFFERENT copy than the one the polar region actually holds --
  //   in which case a minimum-image rule would exclude the wrong one and
  //   silently double-count. Legacy sidesteps the same trap by keying
  //   its own ForegroundTable on explicit (id, na, nb, nc) rather than
  //   on nearest-image. Positions come from JobTopology, which has
  //   already centred them, so an exact-coordinate match is meaningful;
  //   kForegroundMatchTol only absorbs round-off, and is orders of
  //   magnitude below any real inter-segment separation.
  //
  //   Empty by default, which is the plain periodic background.
  EwaldRealSpaceSum(
      const Eigen::Matrix3d& box, const EwaldRegistry& registry, double alpha,
      double thole_a, double r_min, double field_tol,
      double shell_width = 0.945, Index n_max = 15,
      double screening_factor = 6.0,
      const std::vector<std::pair<Index, Eigen::Vector3d>>& foreground = {});

  // Accumulates the total intermolecular real-space field into target's
  // own V()/V_noE() accumulators (via EwaldRealSpaceInteractor, matching
  // its convention), generated by every segment in the registry at charge
  // state source_state, excluding target_segment_id itself (see class
  // documentation on scope). target's position (target.getPos()) is used
  // as the field evaluation point; target itself need not already be
  // registered in registry_. Throws std::runtime_error if convergence is
  // not reached within n_max shells.
  //
  // include_static selects whether each source's PERMANENT multipole
  // field is accumulated alongside its induced-dipole field. It exists
  // for the iterative solver: the permanent contribution does not depend
  // on the induced dipoles, so EwaldPeriodicDipoleOperator's own
  // RawMultiply recomputes an identical constant on every iteration,
  // which its multiply() = RawMultiply(v) - baseline_ then cancels
  // exactly. Passing false there drops that work with no effect on the
  // result (both RawMultiply(v) and baseline_ shift by the same
  // constant). Measured at ~53% of this class's own per-pair real-space
  // cost, i.e. the single largest avoidable cost in the solve.
  //
  // Defaults to true, so the permanent-field pass that builds the
  // solver's own right-hand side (which genuinely needs it) is
  // unaffected.
  template <enum Estatic CE>
  void AddFieldAt(Index target_segment_id, PolarSite& target,
                  EwaldChargeState source_state,
                  bool include_static = true) const;

  // The erfc-screened PERMANENT-multipole interaction energy between one
  // target site and every background source this class would sum a field
  // from -- the same neighbour set, with the same foreground copies
  // suppressed and the same distance cull applied, so the energy and the
  // field are guaranteed to describe the same system.
  //
  // Permanent multipoles on BOTH sides. Neither side's induced dipoles
  // enter: the foreground's are PolarRegion's business (E_polar_ext),
  // and the background's belong to CalcInducedSourceEnergyAt below.
  //
  // Requires the neighbour cache to exist, i.e. AddFieldAt or
  // PrepareNeighborCache must have run for this target first. It does
  // not build the cache itself, because doing so would make an energy
  // query silently expensive and, worse, order-dependent.
  // No target_segment_id parameter, unlike AddFieldAt: the neighbour
  // list is keyed on the target site itself, and the exclusion of the
  // target's own segment is already baked into the cached list, so the
  // id would be dead weight here.
  double CalcStaticEnergyAt(const PolarSite& target,
                            EwaldChargeState source_state) const;

  // The erfc-screened, Thole-damped energy between the BACKGROUND's
  // induced dipoles (sources) and one foreground target's PERMANENT
  // moments -- over exactly the neighbour set CalcStaticEnergyAt uses,
  // for the same reason: the energy and the field must describe the
  // same system.
  //
  // This is the [fg permanent] x [bg induced] corner of the
  // permanent/induced product. See
  // EwaldRealSpaceInteractor::CalcInducedSourceEnergy for why it has no
  // other home, and why it is damped where the permanent energies are
  // not.
  //
  // Same cache requirement as CalcStaticEnergyAt, and for the same
  // reason.
  // POTENTIAL at arbitrary points, batched. Total of both channels, to
  // match EwaldReciprocalSpaceSum::PotentialAtMany and the field
  // AddFieldAt delivers: permanent moments through CalcStaticEnergy,
  // background induced dipoles through CalcInducedSourceEnergy.
  //
  // NO NEIGHBOUR CACHE, deliberately. AddFieldAt keys its cache on the
  // target's ADDRESS and writes an entry on first use, which is right
  // for a fixed set of sites queried repeatedly and wrong for a DFT
  // grid: one entry per point would run to gigabytes, and reusing one
  // probe object across positions would silently hand every later point
  // the first one's neighbour list. Points here are just coordinates, so
  // nothing is keyed and nothing is kept.
  //
  // Its traversal therefore duplicates AddFieldAt's distance cull,
  // foreground suppression and zero-translation self-skip rather than
  // sharing them -- those rules are interleaved with a field-based
  // convergence check there, which has no meaning for a potential.
  // Duplicated rules drift, so this is held to AddFieldAt's own answer
  // by unit_probe_potential_reproduces_the_static_energy and the
  // EwaldRegion case that compares against a direct lattice sum. The
  // shell-convergence early exit is dropped in favour of the geometric
  // cutoff alone, which visits a superset of what the shell search
  // reaches.
  //
  // target_segment_id decides only whether the zero-translation
  // self-pair is skipped, and that skip is disabled for any source with
  // a foreground copy. Pass a foreground segment's id: a point that is
  // not a site of its own wants no self-skip.
  Eigen::VectorXd PotentialAtMany(Index target_segment_id,
                                  const std::vector<Eigen::Vector3d>& points,
                                  EwaldChargeState source_state) const;

  double CalcInducedSourceEnergyAt(const PolarSite& target,
                                   EwaldChargeState source_state) const;

  // Builds the neighbour cache for every target in one serial pass,
  // WITHOUT applying any field (each target's own V()/V_noE() are
  // restored before returning). Exists so callers can parallelize over
  // targets afterwards: neighbor_cache_ and the statistics counters are
  // mutable and written only while a target's list is being built, so
  // once every list exists AddFieldAt is read-only with respect to this
  // object and safe to call concurrently for distinct targets. Calling
  // this is optional -- AddFieldAt still builds its own list on demand
  // -- but a caller that skips it must not run AddFieldAt in parallel.
  void PrepareNeighborCache(
      const std::vector<std::pair<Index, PolarSite*>>& targets,
      EwaldChargeState source_state) const;

  // Neighbour-list statistics over every target whose list has been
  // built so far. entries is the total number of (source segment,
  // translation) pairs kept; culled is how many were rejected by the
  // distance cutoff. Both are zero until the first AddFieldAt call.
  struct NeighborStats {
    Index targets;
    Index entries;
    Index culled;
    Index foreground;
    double entries_per_target() const {
      return targets > 0 ? double(entries) / double(targets) : 0.0;
    }
    double culled_fraction() const {
      const Index seen = entries + culled;
      return seen > 0 ? double(culled) / double(seen) : 0.0;
    }
  };
  NeighborStats GetNeighborStats() const {
    return {cached_targets_, cached_entries_, culled_entries_,
            foreground_entries_};
  }
  double RealSpaceCutoff() const { return real_space_cutoff_; }

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

  // Real-space screening cutoff: the separation beyond which a pair's
  // erfc(alpha*r)-screened contribution is negligible, set to
  // screening_factor / alpha (see the constructor). This is a cost
  // cutoff only: it never extends the sum, and r_min_/field_tol_ still
  // govern how far the shell search goes.
  double real_space_cutoff_;
  // Tolerance for deciding that a shifted source segment coincides with a
  // foreground segment. Round-off only -- see the constructor.
  static constexpr double kForegroundMatchTol = 1e-4;
  // Tolerance for recognising the ZERO lattice translation, used to skip
  // a target's own segment at its own position (intramolecular, and
  // containing the r = 0 self-pair) while keeping its other images. The
  // quantity tested is built from an exactly-zero baseline shift plus an
  // exactly-zero lattice vector, so this only has to absorb round-off;
  // the next-smallest translation is a full lattice vector away.
  static constexpr double kSelfTranslationTol = 1e-8;
  // Foreground copies to suppress, grouped by segment id. Most segments
  // have no entry at all; those that do usually have exactly one.
  std::map<Index, std::vector<Eigen::Vector3d>> foreground_;
  // Neighbour-list statistics, accumulated as the cache is built. The
  // cached list is the real cost driver of the whole solve: every entry
  // is one (source segment, periodic translation) pair, re-evaluated
  // against every one of the target's own sites on every iteration.
  // Reported once per run so the effect of alpha, r_min and the
  // distance cull on that cost is visible directly, rather than being
  // inferred from wall-clock.
  mutable Index cached_targets_ = 0;
  mutable Index cached_entries_ = 0;
  mutable Index culled_entries_ = 0;
  // How many (segment, image) pairs were suppressed as foreground.
  mutable Index foreground_entries_ = 0;
  // Largest site-to-centroid distance over every registered segment,
  // used as the margin when the cutoff (a per-site-pair quantity) is
  // applied at segment granularity.
  double segment_radius_;

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
  // The first tuple element is a direct pointer to the source segment,
  // not its id: resolving an id through registry_ costs a Has() plus a
  // Get(), i.e. two std::map lookups, and this list is walked once per
  // (target site, neighbour) pair -- ~3e7 times per solver iteration on
  // a 1000-segment system, measured at ~59ns per entry against ~3ns for
  // a stored pointer. Legacy's own PolarNb stores the neighbour segment
  // by pointer for the same reason. Safe because the registry outlives
  // this object (see the class documentation) and its segments are
  // stored in a std::map, whose elements keep stable addresses; the same
  // stability argument this cache already relies on for target site
  // addresses being usable as keys.
  mutable std::unordered_map<
      std::pair<const PolarSite*, EwaldChargeState>,
      std::vector<std::tuple<const PolarSegment*, Index, Eigen::Vector3d>>,
      PairHash>
      neighbor_cache_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDREALSPACESUM_H
