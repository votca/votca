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
 * Follows the shell-by-shell scheme of legacy RThread::FP_FieldCalc, but
 * with a simpler convergence criterion: the sum stops once the magnitude
 * of a shell's field contribution drops below an absolute tolerance,
 * provided a minimum radius r_min has been reached (legacy uses a
 * "test dipole" RMS construction instead).
 *
 * Scope: the target's own segment at zero translation (the
 * intramolecular pairs, including the r = 0 self-pair) is excluded; its
 * other periodic images are included. Intramolecular coupling is handled
 * elsewhere, as in the split of eeInteractor between
 * ApplyStaticField/ApplyInducedField and
 * Cholesky_IntraSegment/FillTholeInteraction.
 *
 * Performance: the (source segment, translation) pairs that contribute to
 * a target's converged sum are cached per (target site address, source
 * charge state), so only the first call runs the shell search, as with
 * legacy RThread::PolarNbs(). Without the cache the repeated search made
 * each solver iteration ~8x slower than legacy. The cache assumes that
 * site positions do not change during its lifetime; only induced-dipole
 * values may change.
 *
 * Units: bohr, Hartree atomic units throughout, matching PolarSite and
 * EwaldRealSpaceInteractor.
 */

namespace votca {
namespace xtp {

class EwaldRealSpaceSum {
 public:
  // box: columns are the lattice vectors a, b, c (matching
  //   the Topology::getBox() convention).
  // registry: source of every segment's multipole/induced-dipole state.
  //   Stored by reference; the registry must outlive this object, and
  //   AddFieldAt() reads its current state (e.g. updated induced dipoles).
  // alpha, thole_a: passed through to the underlying
  //   EwaldRealSpaceInteractor.
  // r_min: minimum shell radius (bohr) before convergence may be declared,
  //   regardless of how small any individual shell's contribution is.
  // field_tol: convergence threshold (atomic units) on a shell's
  //   contribution to the field magnitude.
  // shell_width: radial bin width (bohr) used to group periodic images
  //   into shells.
  // n_max: safety cap on the lattice-translation search box (translations
  //   with |na|,|nb|,|nc| <= n_max are considered candidates). If
  //   convergence is not reached within this box, AddFieldAt() throws.
  // screening_factor: sets the real-space distance cutoff to
  //   screening_factor / alpha. Beyond it a pair contributes at most
  //   ~erfc(screening_factor); the default 6.0 gives ~2e-17 and discards
  //   78.7% of the pairs on a 1000-segment box (4.5x faster real-space
  //   sum, converged dipoles unchanged in all printed digits). The number
  //   of retained pairs grows as the cube of the cutoff. Values below ~4
  //   affect results at the 1e-11 level.
  //   Kept as the last parameter so that existing positional calls are
  //   not silently reinterpreted; add new parameters at the end.
  //
  // foreground: (segment id, centroid position) pairs treated explicitly
  //   by a polar region and therefore removed from this sum (the MM/MM or
  //   QM/MM carve-out). Each entry suppresses only the periodic copy at
  //   the given position; all other images remain in the background.
  //
  //   Copies are matched by position, not by minimum image: two
  //   foreground segments can be up to twice the region cutoff apart, so
  //   the minimum image may differ from the copy held by the polar region
  //   (legacy likewise keys its ForegroundTable on (id, na, nb, nc)).
  //   Positions come from JobTopology and are already centred, so
  //   kForegroundMatchTol only absorbs round-off.
  //
  //   Empty by default, which is the plain periodic background.
  EwaldRealSpaceSum(
      const Eigen::Matrix3d& box, const EwaldRegistry& registry, double alpha,
      double thole_a, double r_min, double field_tol,
      double shell_width = 0.945, Index n_max = 15,
      double screening_factor = 6.0,
      const std::vector<std::pair<Index, Eigen::Vector3d>>& foreground = {});

  // Accumulates the total intermolecular real-space field into target's
  // V()/V_noE() accumulators (EwaldRealSpaceInteractor convention) from
  // every segment registered at source_state, excluding the
  // zero-translation copy of target_segment_id (see class documentation).
  // The field is evaluated at target.getPos(); target need not be
  // registered. Throws std::runtime_error if convergence is not reached
  // within n_max shells.
  //
  // include_static = false skips the permanent-multipole field and
  // accumulates only the induced-dipole field. EwaldPeriodicDipoleOperator
  // uses this, since the permanent field is a constant that cancels in
  // RawMultiply(v) - baseline_; it is ~53% of the per-pair cost.
  template <enum Estatic CE>
  void AddFieldAt(Index target_segment_id, PolarSite& target,
                  EwaldChargeState source_state,
                  bool include_static = true) const;

  // AddFieldAt for a whole batch, in parallel. Every missing neighbour
  // list is built first by PrepareNeighborCache (parallel search, serial
  // insertion), after which the field pass only reads this object.
  // Targets are spread over threads by site: repeated entries of the same
  // site are applied by one thread, in input order. The result is
  // bit-identical to calling AddFieldAt on every entry in turn: each
  // target sees the same (source, translation) list, applied in the same
  // order.
  //
  // Use this (or PrepareNeighborCache first) rather than an OpenMP loop
  // around AddFieldAt: AddFieldAt inserts into neighbor_cache_ on a miss,
  // and concurrent insertion into an unordered_map is a data race.
  template <enum Estatic CE>
  void AddFieldAtMany(const std::vector<std::pair<Index, PolarSite*>>& targets,
                      EwaldChargeState source_state,
                      bool include_static = true) const;

  // erfc-screened interaction energy between the PERMANENT multipoles of
  // one target site and of every background source, over the same cached
  // neighbour list (foreground suppression and distance cull included)
  // that AddFieldAt uses, so energy and field describe the same system.
  // Induced dipoles do not enter: the foreground's are handled by
  // PolarRegion (E_polar_ext), the background's by
  // CalcInducedSourceEnergyAt.
  //
  // Requires the neighbour list of this target (AddFieldAt or
  // PrepareNeighborCache must have run); it is not built here, to keep
  // the energy query cheap and independent of call order. The list
  // already encodes the own-segment exclusion, so no segment id is
  // needed.
  double CalcStaticEnergyAt(const PolarSite& target,
                            EwaldChargeState source_state) const;

  // POTENTIAL at arbitrary points, batched. Total of both channels, to
  // match EwaldReciprocalSpaceSum::PotentialAtMany and the field
  // AddFieldAt delivers: permanent moments through CalcStaticEnergy,
  // background induced dipoles through CalcInducedSourceEnergy.
  //
  // Uses no neighbour cache: the cache is keyed on site addresses, and
  // one entry per grid point (e.g. a DFT grid) would be prohibitively
  // large. The distance cull, foreground suppression and zero-translation
  // self-skip of AddFieldAt are therefore duplicated here; consistency is
  // tested by unit_probe_potential_reproduces_the_static_energy and the
  // EwaldRegion comparison with a direct lattice sum. Instead of the
  // shell-convergence exit only the geometric cutoff is used, which
  // visits a superset of the shell search.
  //
  // target_segment_id decides only whether the zero-translation
  // self-pair is skipped, and that skip is disabled for any source with
  // a foreground copy. Pass a foreground segment's id: a point that is
  // not a site of its own wants no self-skip.
  Eigen::VectorXd PotentialAtMany(Index target_segment_id,
                                  const std::vector<Eigen::Vector3d>& points,
                                  EwaldChargeState source_state) const;

  // erfc-screened, Thole-damped energy between the BACKGROUND induced
  // dipoles (sources) and the PERMANENT moments of one foreground target,
  // over the same neighbour list as CalcStaticEnergyAt (same cache
  // requirement). This is the [fg permanent] x [bg induced] term; see
  // EwaldRealSpaceInteractor::CalcInducedSourceEnergy for why it is
  // damped where the permanent energies are not.
  double CalcInducedSourceEnergyAt(const PolarSite& target,
                                   EwaldChargeState source_state) const;

  // Builds the neighbour list of every target that does not have one
  // yet, WITHOUT applying any field: the search runs on a scratch copy of
  // each target, so the targets themselves are not touched. Targets that
  // already have a list (and repeated entries) are skipped.
  //
  // The searches run in parallel; the lists and statistics are then
  // inserted serially. Afterwards AddFieldAt does not modify this object
  // and may be called concurrently for distinct targets. Without this
  // call AddFieldAt builds lists on demand and must not run in parallel.
  //
  // converge_static_field selects the field whose shell contribution
  // terminates the search, as include_static does for AddFieldAt: true
  // for permanent plus induced, false for induced only. Pass the value
  // the first AddFieldAt call would use, so that the list is identical
  // to the on-demand one.
  void PrepareNeighborCache(
      const std::vector<std::pair<Index, PolarSite*>>& targets,
      EwaldChargeState source_state, bool converge_static_field = false) const;

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

  using NeighborList =
      std::vector<std::tuple<const PolarSegment*, Index, Eigen::Vector3d>>;

  // Per-search statistics, accumulated locally so that a search can run
  // without writing to this object, and added to the member counters by
  // whoever inserts the list.
  struct SearchCounts {
    Index culled = 0;
    Index foreground = 0;
  };

  // The shell-by-shell convergence search for one target. Accumulates
  // the field into scratch (a copy of the target: the convergence test
  // needs a field to watch) and returns every (source, translation) pair
  // it applied, in the order it applied them. Writes nothing in this
  // object, so concurrent calls are safe.
  template <enum Estatic CE>
  NeighborList SearchNeighbors(Index target_segment_id, PolarSite& scratch,
                               EwaldChargeState source_state,
                               bool include_static, SearchCounts& counts) const;

  // Applies a neighbour list to target: the fast path of AddFieldAt.
  template <enum Estatic CE>
  void ApplyNeighborList(const NeighborList& list, PolarSite& target,
                         bool include_static) const;

  // PrepareNeighborCache for the search AddFieldAt<CE>(..., include_static)
  // would run on a miss, so a list built here is the one the on-demand
  // path would have built. AddFieldAtMany<CE> uses this directly.
  template <enum Estatic CE>
  void PrepareNeighborCacheImpl(
      const std::vector<std::pair<Index, PolarSite*>>& targets,
      EwaldChargeState source_state, bool include_static) const;

  // Inserts one list and its statistics. Not thread-safe.
  void InsertNeighborList(
      const std::pair<const PolarSite*, EwaldChargeState>& key,
      NeighborList list, const SearchCounts& counts) const;

  // std::pair has no std::hash specialization; combine the hashes of
  // the pointer and the enum (hashable since C++14).
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
  // number of entries determines the cost per solver iteration; it is
  // reported so that the effect of alpha, r_min and the cull is visible.
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
  // target position), reused by every AddFieldAt() call.
  std::vector<Translation> translations_;

  // Neighbour list (source segment, translation) of each target's
  // converged sum, analogous to legacy RThread::PolarNbs(). Keyed by
  // (target site address, source charge state), since the set of
  // registered sources depends on the charge state.
  // Sources are stored by pointer rather than id: an id lookup costs two
  // std::map lookups (~59 ns vs ~3 ns per entry, ~3e7 entries per
  // iteration for 1000 segments). This is safe because the registry
  // outlives this object and stores segments in a std::map, whose
  // elements have stable addresses.
  mutable std::unordered_map<std::pair<const PolarSite*, EwaldChargeState>,
                             NeighborList, PairHash>
      neighbor_cache_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDREALSPACESUM_H
