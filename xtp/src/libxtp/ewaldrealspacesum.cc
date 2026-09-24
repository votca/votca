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

// Standard includes
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Local VOTCA includes
#include "votca/xtp/ewaldrealspacesum.h"

namespace votca {
namespace xtp {

namespace {
// BUG FIX (this session): PolarSegment::getPos() (via AtomContainer<T>::
// calcPos()) is a MASS-weighted center of mass. Legacy's own equivalent
// (PolarSeg::CalcPos()) is a plain, UNWEIGHTED arithmetic mean of site
// positions instead -- confirmed by direct comparison of the two
// implementations. Since this position feeds directly into the
// minimum-image PBC wrap below (raw_offset = target - source's own
// representative position), even a small difference between the two
// definitions can flip which periodic image is chosen as the true
// minimum for a given pair -- a discrete, not continuous, effect,
// consistent with the scattered (not uniformly scaled) per-site
// mismatch pattern that motivated this fix, rather than a smooth
// site-by-site drift.
//
// Deliberately NOT fixed by changing calcPos() itself: that's a shared
// AtomContainer<T> method used well beyond this file (md2qmengine,
// segmentmapper), where a genuine physical mass-weighted center of mass
// may be exactly what's wanted. This local helper instead gives THIS
// file its own, legacy-matching definition, without touching shared
// code any other caller depends on.
Eigen::Vector3d UnweightedCentroid(const PolarSegment& seg) {
  Eigen::Vector3d pos = Eigen::Vector3d::Zero();
  Index n = 0;
  for (const PolarSite& site : seg) {
    pos += site.getPos();
    ++n;
  }
  if (n > 0) {
    pos /= double(n);
  }
  return pos;
}
}  // namespace

EwaldRealSpaceSum::EwaldRealSpaceSum(const Eigen::Matrix3d& box,
                                     const EwaldRegistry& registry,
                                     double alpha, double thole_a,
                                     double r_min, double field_tol,
                                     double shell_width, Index n_max,
                                     double screening_factor,
                                     const std::vector<
                                         std::pair<Index, Eigen::Vector3d>>&
                                         foreground)
    : real_space_cutoff_(screening_factor / alpha),
      segment_radius_(0.0),
      box_(box),
      registry_(registry),
      interactor_(alpha, thole_a),
      r_min_(r_min),
      field_tol_(field_tol),
      shell_width_(shell_width),
      n_max_(n_max) {
  translations_ = GenerateSortedTranslations();

  for (const auto& entry : foreground) {
    foreground_[entry.first].push_back(entry.second);
  }

  // Largest site-to-centroid distance anywhere in the registry -- the
  // margin the segment-granular distance cull needs (see AddFieldAt).
  // Computed once here rather than per call: the registry's geometry is
  // fixed for this object's lifetime, the same assumption the neighbor
  // cache already relies on.
  for (Index id : registry_.AllIds()) {
    for (EwaldChargeState state :
         {EwaldChargeState::Neutral, EwaldChargeState::Electron,
          EwaldChargeState::Hole}) {
      if (!registry_.Has(id, state)) {
        continue;
      }
      const PolarSegment& segment = registry_.Get(id, state);
      const Eigen::Vector3d centroid = UnweightedCentroid(segment);
      for (const PolarSite& site : segment) {
        segment_radius_ =
            std::max(segment_radius_, (site.getPos() - centroid).norm());
      }
    }
  }
}

std::vector<EwaldRealSpaceSum::Translation>
EwaldRealSpaceSum::GenerateSortedTranslations() const {
  std::vector<Translation> translations;
  const Eigen::Vector3d a = box_.col(0);
  const Eigen::Vector3d b = box_.col(1);
  const Eigen::Vector3d c = box_.col(2);

  for (Index na = -n_max_; na <= n_max_; ++na) {
    for (Index nb = -n_max_; nb <= n_max_; ++nb) {
      for (Index nc = -n_max_; nc <= n_max_; ++nc) {
        Eigen::Vector3d t = double(na) * a + double(nb) * b + double(nc) * c;
        translations.push_back({t, t.norm()});
      }
    }
  }

  std::sort(translations.begin(), translations.end(),
           [](const Translation& x, const Translation& y) {
             return x.r < y.r;
           });
  return translations;
}

template <enum Estatic CE>
void EwaldRealSpaceSum::AddFieldAt(Index target_segment_id, PolarSite& target,
                                   EwaldChargeState source_state,
                                   bool include_static) const {
  const std::pair<const PolarSite*, EwaldChargeState> cache_key(&target,
                                                                source_state);
  auto cached = neighbor_cache_.find(cache_key);
  if (cached != neighbor_cache_.end()) {
    // Fast path: the real, geometry-determined neighbor set for this
    // target was already found by an earlier call (see this class's own
    // documentation for why reusing it across calls for the same target
    // is safe). Just apply every cached (source, translation) pair
    // directly -- no shell-by-shell search, no scan of
    // registry_.AllIds() members that turned out not to matter.
    for (const auto& entry : cached->second) {
      // Stored by pointer rather than by id -- see neighbor_cache_'s own
      // declaration for why. The Has()/Get() guard the id-based version
      // needed is gone with it: the cache is keyed on source_state, so
      // every entry in THIS list was registered at THIS state when the
      // list was built, and segment addresses are stable for the
      // registry's lifetime.
      const PolarSegment& source_segment = *std::get<0>(entry);
      const Index translation_idx = std::get<1>(entry);
      const Eigen::Vector3d& baseline_shift = std::get<2>(entry);
      const Eigen::Vector3d t =
          baseline_shift + translations_[translation_idx].t;
      for (const PolarSite& source_site : source_segment) {
        // Shift passed through rather than applied to a copy of the
        // source site -- see ApplyStaticField's own source_shift note.
        if (include_static) {
          interactor_.ApplyStaticField<PolarSite, CE>(source_site, target, t);
        }
        interactor_.ApplyInducedField<CE>(source_site, target, t);
      }
    }
    return;
  }

  // Slow path: first call for this target. Runs the full convergence
  // search, exactly as before, but additionally records every (source,
  // translation) pair actually visited, so every subsequent call for
  // this same target can skip straight to the fast path above.
  std::vector<std::tuple<const PolarSegment*, Index, Eigen::Vector3d>>
      visited_pairs;

  // See the cull inside the shell loop below for what this is and why.
  const double cutoff_with_margin =
      real_space_cutoff_ + 2.0 * segment_radius_;

  Index shell_start = 0;
  double shell_edge = 0.0;
  bool converged = false;

  while (shell_start < Index(translations_.size())) {
    // Advance shell_edge by shell_width_ until it covers at least the next
    // not-yet-processed translation, then collect every translation up to
    // that edge into this shell. Since translations_ is sorted by
    // distance, this always yields a contiguous, radially-ordered shell.
    shell_edge = translations_[shell_start].r +
                (shell_edge > translations_[shell_start].r ? 0.0 : shell_width_);
    Index shell_end = shell_start;
    while (shell_end < Index(translations_.size()) &&
          translations_[shell_end].r <= shell_edge) {
      ++shell_end;
    }

    // Field is written into target.V()/V_noE() by ApplyStaticField/
    // ApplyInducedField as a side effect; recover this shell's own
    // contribution by differencing target's accumulator before and after
    // processing the whole shell.
    const Eigen::Vector3d before_shell =
        (CE == Estatic::noE_V) ? target.V_noE() : target.V();

    for (Index source_id : registry_.AllIds()) {
      // No per-segment skip: a target's own segment is treated like any
      // other. Its COINCIDENT copy is dropped by the foreground
      // suppression below; its other lattice images are ordinary
      // background molecules and stay. Legacy agrees -- SetupMidground
      // excludes by (segment id, na, nb, nc), not by segment id.
      if (!registry_.Has(source_id, source_state)) {
        continue;
      }
      const PolarSegment& source_segment =
          registry_.Get(source_id, source_state);

      // Per-source baseline wrap (computed once per source_id, not per
      // shell-translation index, since it depends only on the target
      // and source_segment's own representative positions, not on
      // which small shell-translation t is currently being applied).
      //
      // BUG FIX (this session, found via a real production comparison
      // against legacy): translations_ is sorted and shell-converged
      // by |t| alone (the raw lattice-translation magnitude), but the
      // quantity that actually determines a pair's true separation is
      // |raw_offset + t|, not |t| by itself. A translation with a
      // LARGE |t| can still be the one that gives the SMALLEST actual
      // distance, if it happens to cancel a large raw target-source
      // offset -- and such a translation is exactly what the
      // shell-convergence criterion (stop once shell radius >= r_min_)
      // can miss entirely, since it never even considers translations
      // whose OWN magnitude exceeds r_min_, regardless of what they'd
      // produce once combined with the raw offset. Confirmed on a real
      // pair: legacy's own minimum-image search (round-based,
      // independent of any shell-magnitude cutoff) found r=13.6 bohr
      // via t=(0,-1,-1)*box (|t|~97 bohr, comfortably past r_min_ in
      // this run's own ~75.6 bohr setting) for a pair this class's own
      // unfixed search reported no closer than r=55.4 bohr on, because
      // it never got there.
      //
      // Fixed by mirroring legacy's own two-step structure exactly:
      // wrap the raw (large, unbounded) target-source separation into
      // its true minimum image FIRST (fractional-coordinate rounding,
      // correct for any box_ shape, not just orthorhombic), THEN run
      // the existing small, bounded shell search of translations_ on
      // top of that already-small baseline -- exactly as legacy's own
      // dr12_pbc + L structure does. The shell search only ever needs
      // to explore the neighborhood right around the true minimum
      // image, not the raw, unbounded separation, so its own
      // shell_radius >= r_min_ stopping criterion is now sound (it was
      // never wrong in isolation -- it was being applied to the wrong
      // starting point).
      const Eigen::Vector3d source_centroid =
          UnweightedCentroid(source_segment);
      const Eigen::Vector3d raw_offset = target.getPos() - source_centroid;
      const Eigen::Vector3d frac = box_.inverse() * raw_offset;
      const Eigen::Vector3d wrapped_frac =
          frac - frac.array().round().matrix();
      const Eigen::Vector3d min_image_offset = box_ * wrapped_frac;
      const Eigen::Vector3d baseline_shift =
          raw_offset - min_image_offset;

      for (Index idx = shell_start; idx < shell_end; ++idx) {
        // Distance cull. The erfc(alpha*r) screening means a pair's
        // contribution falls off as fast as erfc does: at alpha*r = 6 a
        // single pair contributes ~1e-20 of the total field, and the
        // whole remaining tail sums to far below field_tol_. Without
        // this test every registered segment is evaluated at every
        // accepted translation, however far away it is -- measured on a
        // real 1000-segment system as 14985 (segment, translation)
        // entries per target, of which only ~22% lie inside the cutoff.
        // That factor of ~4.6, not any per-pair cost, was the whole
        // remaining real-space gap against legacy (whose own
        // PolarNbs list is distance-built, so it never had the excess).
        //
        // r_min_ does NOT already do this: it bounds the shell search by
        // TRANSLATION magnitude |t|, which is unrelated to how far a
        // given source segment is once that translation is applied.
        //
        // The true separation for this (segment, translation) pair is
        // |min_image_offset - translations_[idx].t| -- the same quantity
        // the interactor will form, since it evaluates
        // target - (source + t) with t = baseline_shift +
        // translations_[idx].t and baseline_shift = raw_offset -
        // min_image_offset. Compared at segment granularity, so
        // segment_radius_ (the largest site-to-centroid distance in the
        // registry) is added twice as a margin: no individual site pair
        // can then be closer than the cutoff while its segment pair is
        // culled.
        const double pair_distance =
            (min_image_offset - translations_[idx].t).norm();
        if (pair_distance > cutoff_with_margin) {
          ++culled_entries_;
          continue;
        }
        const Eigen::Vector3d t = baseline_shift + translations_[idx].t;

        // Foreground suppression. This one periodic copy of this segment
        // is handled explicitly elsewhere (a polar region), so it must
        // not also appear in the periodic background -- see the
        // constructor's own foreground documentation. Only the copy
        // sitting at the recorded position is dropped; the segment's
        // other lattice images stay in the sum.
        if (!foreground_.empty()) {
          auto fg = foreground_.find(source_id);
          if (fg != foreground_.end()) {
            const Eigen::Vector3d shifted_centroid = source_centroid + t;
            bool suppressed = false;
            for (const Eigen::Vector3d& fg_pos : fg->second) {
              if ((shifted_centroid - fg_pos).norm() < kForegroundMatchTol) {
                suppressed = true;
                break;
              }
            }
            if (suppressed) {
              ++foreground_entries_;
              continue;
            }
          }
        }

        // The target's own segment at ZERO translation is the segment
        // itself: intramolecular pairs plus the r = 0 self-pair, never
        // part of this sum (AddIntraSegmentCoupling owns the induced
        // side, the permanent side has its own compensation pass).
        //
        // Only where the segment has no recorded foreground copy. Where
        // it does, the suppression above already decides which copy is
        // carved out, and that need not be the t = 0 one: a foreground
        // segment sitting at a nonzero image has its t = 0 copy as a
        // genuine neighbour. Omitting this guard made the field infinite
        // in every background solve.
        const bool source_has_foreground =
            !foreground_.empty() && foreground_.count(source_id) > 0;
        if (!source_has_foreground && source_id == target_segment_id &&
            t.squaredNorm() < kSelfTranslationTol * kSelfTranslationTol) {
          continue;
        }

        visited_pairs.emplace_back(&source_segment, idx, baseline_shift);
        for (const PolarSite& source_site : source_segment) {
          if (include_static) {
            interactor_.ApplyStaticField<PolarSite, CE>(source_site, target, t);
          }
          interactor_.ApplyInducedField<CE>(source_site, target, t);
        }
      }
    }

    const Eigen::Vector3d after_shell =
        (CE == Estatic::noE_V) ? target.V_noE() : target.V();
    const Eigen::Vector3d shell_field = after_shell - before_shell;

    const double shell_radius =
        translations_[shell_end > shell_start ? shell_end - 1 : shell_start]
            .r;
    if (shell_radius >= r_min_ && shell_field.norm() < field_tol_) {
      converged = true;
      break;
    }

    shell_start = shell_end;
  }

  if (!converged) {
    std::stringstream message;
    message << "EwaldRealSpaceSum: real-space sum for segment "
            << target_segment_id
            << " did not converge within the n_max lattice-translation "
               "search box; increase n_max or field_tol.";
    throw std::runtime_error(message.str());
  }

  ++cached_targets_;
  cached_entries_ += Index(visited_pairs.size());
  neighbor_cache_.emplace(cache_key, std::move(visited_pairs));
}

double EwaldRealSpaceSum::CalcStaticEnergyAt(
    const PolarSite& target, EwaldChargeState source_state) const {
  // See this method's own declaration for what is and is not included.
  const std::pair<const PolarSite*, EwaldChargeState> cache_key(&target,
                                                                 source_state);
  auto cached = neighbor_cache_.find(cache_key);
  if (cached == neighbor_cache_.end()) {
    throw std::runtime_error(
        "EwaldRealSpaceSum::CalcStaticEnergyAt: no neighbour list for this "
        "target. AddFieldAt or PrepareNeighborCache must run first -- this "
        "method deliberately does not build one, so that an energy query "
        "cannot silently become the expensive shell search.");
  }

  double energy = 0.0;
  for (const auto& entry : cached->second) {
    const PolarSegment& source_segment = *std::get<0>(entry);
    const Index translation_idx = std::get<1>(entry);
    const Eigen::Vector3d& baseline_shift = std::get<2>(entry);
    const Eigen::Vector3d t =
        baseline_shift + translations_[translation_idx].t;
    for (const PolarSite& source_site : source_segment) {
      energy += interactor_.CalcStaticEnergy<PolarSite, PolarSite>(
          source_site, target, t);
    }
  }
  return energy;
}

Eigen::VectorXd EwaldRealSpaceSum::PotentialAtMany(
    Index target_segment_id, const std::vector<Eigen::Vector3d>& points,
    EwaldChargeState source_state) const {
  // Sources flattened once: AllIds() returns by value and the centroids
  // are fixed, so recomputing either per point would cost more than the
  // sum itself on a grid.
  std::vector<const PolarSegment*> segments;
  std::vector<Index> segment_ids;
  std::vector<Eigen::Vector3d> centroids;
  for (Index source_id : registry_.AllIds()) {
    if (!registry_.Has(source_id, source_state)) {
      continue;
    }
    const PolarSegment& segment = registry_.Get(source_id, source_state);
    segments.push_back(&segment);
    segment_ids.push_back(source_id);
    centroids.push_back(UnweightedCentroid(segment));
  }

  const double cutoff_with_margin = real_space_cutoff_ + 2.0 * segment_radius_;
  const Index n_points = Index(points.size());
  const Index n_sources = Index(segments.size());
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(n_points);

  // Hoisted out of both loops below: the box is fixed for this object's
  // lifetime, and inverting it per (point, source) is n_points * n_sources
  // 3x3 inversions for one matrix.
  const Eigen::Matrix3d box_inv = box_.inverse();

#pragma omp parallel for schedule(dynamic, 8)
  for (Index p = 0; p < n_points; ++p) {
    const Eigen::Vector3d& point = points[std::size_t(p)];

    // Unit test charge: every energy routine here reduces to
    // q*phi - mu.E, so with q = 1 and mu = 0 what comes back is phi.
    // Built per point on the stack, which is safe only because nothing
    // below keys anything on its address.
    PolarSite probe(0, "H", point);
    probe.Reset();
    probe.setCharge(1.0);
    probe.setStaticDipole(Eigen::Vector3d::Zero());
    probe.setInduced_Dipole(Eigen::Vector3d::Zero());

    double acc = 0.0;
    for (Index s_i = 0; s_i < n_sources; ++s_i) {
      const PolarSegment& source_segment = *segments[std::size_t(s_i)];
      const Index source_id = segment_ids[std::size_t(s_i)];
      const Eigen::Vector3d& source_centroid = centroids[std::size_t(s_i)];

      // Minimum image first, shell search on top of it -- see AddFieldAt's
      // own account of why the raw separation cannot be handed to a
      // search bounded by |t|.
      const Eigen::Vector3d raw_offset = point - source_centroid;
      const Eigen::Vector3d frac = box_inv * raw_offset;
      const Eigen::Vector3d wrapped_frac =
          frac - frac.array().round().matrix();
      const Eigen::Vector3d min_image_offset = box_ * wrapped_frac;
      const Eigen::Vector3d baseline_shift = raw_offset - min_image_offset;

      const auto fg = foreground_.find(source_id);
      const bool source_has_foreground = fg != foreground_.end();

      // STOP, don't skip. translations_ is sorted by |t| ascending, and
      // |min_image_offset - t| >= |t| - |min_image_offset|, so once |t|
      // passes this limit every remaining translation is culled as well.
      // The bound is exact, not a heuristic: nothing inside the cutoff can
      // sit beyond it.
      //
      // This matters far more here than in the site-based paths above,
      // which never walk the whole list -- their shell-convergence search
      // hands them a [shell_start, shell_end) window. Without the break
      // this loop runs (2*n_max+1)^3 translations, 29791 at the default
      // n_max of 15, for EVERY (point, source) pair. On a DFT integration
      // grid against a 1000-segment registry that is of order 1e12
      // distance evaluations and hours of wall time; with it the surviving
      // range is a couple of lattice shells.
      const double t_limit = cutoff_with_margin + min_image_offset.norm();

      for (std::size_t idx = 0; idx < translations_.size(); ++idx) {
        if (translations_[idx].r > t_limit) {
          break;
        }
        const double pair_distance =
            (min_image_offset - translations_[idx].t).norm();
        if (pair_distance > cutoff_with_margin) {
          continue;
        }
        const Eigen::Vector3d t = baseline_shift + translations_[idx].t;

        // Foreground suppression: this one periodic copy is handled
        // explicitly elsewhere, so it must not also appear here.
        if (source_has_foreground) {
          const Eigen::Vector3d shifted_centroid = source_centroid + t;
          bool suppressed = false;
          for (const Eigen::Vector3d& fg_pos : fg->second) {
            if ((shifted_centroid - fg_pos).norm() < kForegroundMatchTol) {
              suppressed = true;
              break;
            }
          }
          if (suppressed) {
            continue;
          }
        }

        if (!source_has_foreground && source_id == target_segment_id &&
            t.squaredNorm() < kSelfTranslationTol * kSelfTranslationTol) {
          continue;
        }

        for (const PolarSite& source_site : source_segment) {
          acc += interactor_.CalcStaticEnergy<PolarSite, PolarSite>(
              source_site, probe, t);
          acc += interactor_.CalcInducedSourceEnergy(source_site, probe, t);
        }
      }
    }
    phi[p] = acc;
  }
  return phi;
}

double EwaldRealSpaceSum::CalcInducedSourceEnergyAt(
    const PolarSite& target, EwaldChargeState source_state) const {
  // Deliberately a near-copy of CalcStaticEnergyAt above rather than a
  // shared template over the interactor call: the two differ only in
  // which moments they contract, but they are validated separately and
  // against different legacy channels (_pp and half of _pu), and a
  // shared body would let a change to one silently move the other.
  const std::pair<const PolarSite*, EwaldChargeState> cache_key(&target,
                                                                 source_state);
  auto cached = neighbor_cache_.find(cache_key);
  if (cached == neighbor_cache_.end()) {
    throw std::runtime_error(
        "EwaldRealSpaceSum::CalcInducedSourceEnergyAt: no neighbour list for "
        "this target. AddFieldAt or PrepareNeighborCache must run first -- "
        "this method deliberately does not build one, so that an energy "
        "query cannot silently become the expensive shell search.");
  }

  double energy = 0.0;
  for (const auto& entry : cached->second) {
    const PolarSegment& source_segment = *std::get<0>(entry);
    const Index translation_idx = std::get<1>(entry);
    const Eigen::Vector3d& baseline_shift = std::get<2>(entry);
    const Eigen::Vector3d t =
        baseline_shift + translations_[translation_idx].t;
    for (const PolarSite& source_site : source_segment) {
      energy += interactor_.CalcInducedSourceEnergy(source_site, target, t);
    }
  }
  return energy;
}

void EwaldRealSpaceSum::PrepareNeighborCache(
    const std::vector<std::pair<Index, PolarSite*>>& targets,
    EwaldChargeState source_state) const {
  // See this method's own declaration for what this is for. The search
  // and the field accumulation share one code path in AddFieldAt, so
  // the list is built by running it and then undoing its effect on the
  // target, rather than by duplicating the shell-search logic here
  // where the two copies could drift apart.
  for (const auto& entry : targets) {
    PolarSite& target = *entry.second;
    const Eigen::Vector3d saved_V = target.V();
    const Eigen::Vector3d saved_V_noE = target.V_noE();
    AddFieldAt<Estatic::V>(entry.first, target, source_state, false);
    target.V() = saved_V;
    target.V_noE() = saved_V_noE;
  }
}

template void EwaldRealSpaceSum::AddFieldAt<Estatic::V>(
    Index, PolarSite&, EwaldChargeState, bool) const;
template void EwaldRealSpaceSum::AddFieldAt<Estatic::noE_V>(
    Index, PolarSite&, EwaldChargeState, bool) const;




}  // namespace xtp
}  // namespace votca
