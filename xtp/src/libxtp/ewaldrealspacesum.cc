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

EwaldRealSpaceSum::EwaldRealSpaceSum(const Eigen::Matrix3d& box,
                                     const EwaldRegistry& registry,
                                     double alpha, double thole_a,
                                     double r_min, double field_tol,
                                     double shell_width, Index n_max)
    : box_(box),
      registry_(registry),
      interactor_(alpha, thole_a),
      r_min_(r_min),
      field_tol_(field_tol),
      shell_width_(shell_width),
      n_max_(n_max) {
  translations_ = GenerateSortedTranslations();
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
                                   EwaldChargeState source_state) const {
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
      const Index source_id = std::get<0>(entry);
      const Index translation_idx = std::get<1>(entry);
      const Eigen::Vector3d& baseline_shift = std::get<2>(entry);
      if (!registry_.Has(source_id, source_state)) {
        // A source segment present at cache-build time might not be
        // registered at this call's own source_state (e.g. built while
        // iterating a different EwaldChargeState) -- skip rather than
        // let Get() throw, matching the original search's own
        // !registry_.Has(...) guard.
        continue;
      }
      const PolarSegment& source_segment =
          registry_.Get(source_id, source_state);
      const Eigen::Vector3d t =
          baseline_shift + translations_[translation_idx].t;
      for (const PolarSite& source_site : source_segment) {
        PolarSite shifted = source_site;
        shifted.setPos(source_site.getPos() + t);

        interactor_.ApplyStaticField<PolarSite, CE>(shifted, target);
        interactor_.ApplyInducedField<CE>(shifted, target);
      }
    }
    return;
  }

  // Slow path: first call for this target. Runs the full convergence
  // search, exactly as before, but additionally records every (source,
  // translation) pair actually visited, so every subsequent call for
  // this same target can skip straight to the fast path above.
  std::vector<std::tuple<Index, Index, Eigen::Vector3d>> visited_pairs;

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
      if (source_id == target_segment_id) {
        // Intermolecular scope only -- see class documentation. This
        // excludes every site of the target's own segment, at every
        // translation (including nonzero ones, i.e. that segment's own
        // periodic images), not just the zero-translation case.
        continue;
      }
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
      const Eigen::Vector3d raw_offset =
          target.getPos() - source_segment.getPos();
      const Eigen::Vector3d frac = box_.inverse() * raw_offset;
      const Eigen::Vector3d wrapped_frac =
          frac - frac.array().round().matrix();
      const Eigen::Vector3d min_image_offset = box_ * wrapped_frac;
      const Eigen::Vector3d baseline_shift =
          raw_offset - min_image_offset;

      for (Index idx = shell_start; idx < shell_end; ++idx) {
        visited_pairs.emplace_back(source_id, idx, baseline_shift);
        const Eigen::Vector3d t = baseline_shift + translations_[idx].t;
        for (const PolarSite& source_site : source_segment) {
          PolarSite shifted = source_site;
          shifted.setPos(source_site.getPos() + t);

          interactor_.ApplyStaticField<PolarSite, CE>(shifted, target);
          interactor_.ApplyInducedField<CE>(shifted, target);
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

  neighbor_cache_.emplace(cache_key, std::move(visited_pairs));
}

template void EwaldRealSpaceSum::AddFieldAt<Estatic::V>(
    Index, PolarSite&, EwaldChargeState) const;
template void EwaldRealSpaceSum::AddFieldAt<Estatic::noE_V>(
    Index, PolarSite&, EwaldChargeState) const;

void EwaldRealSpaceSum::DumpNeighborListAppend(
    Index target_segment_id, const PolarSite& target,
    EwaldChargeState source_state, const std::string& filename) const {
  // See this method's own declaration for what this is for. Mirrors
  // AddFieldAt's own slow-path shell search exactly (same
  // shell_start/shell_edge/shell_end advancing logic, same
  // registry_.AllIds() scan, same target_segment_id exclusion) --
  // deliberately not reusing or populating neighbor_cache_, so this
  // dump reflects a genuinely independent re-run of the search, not a
  // cached result from an earlier AddFieldAt call for this same
  // target. Always appends -- see this method's own declaration for
  // why.
  //
  // BUG FIX (this session, after a real production run showed this
  // giving r=0 for every single source, universally): an earlier
  // version of this dumped translations_[idx].r and .t directly --
  // the raw lattice-vector magnitude, which trivially includes
  // t=(0,0,0) (r=0) for every source regardless of the true
  // target-source separation, since it never added the actual
  // (source position + t) - target position offset the way
  // AddFieldAt's own field computation does. Fixed here to compute
  // the real distance for each translation, using the source
  // segment's own first site as its representative position --
  // matching this dump's own segment-level granularity (mirroring
  // legacy's own PolarNb, which likewise wraps a whole segment, not
  // individual sites).
  std::ofstream dump(filename, std::ios::app);
  dump.precision(15);

  Index shell_start = 0;
  double shell_edge = 0.0;
  bool converged = false;

  while (shell_start < Index(translations_.size())) {
    shell_edge = translations_[shell_start].r +
                (shell_edge > translations_[shell_start].r ? 0.0
                                                            : shell_width_);
    Index shell_end = shell_start;
    while (shell_end < Index(translations_.size()) &&
          translations_[shell_end].r <= shell_edge) {
      ++shell_end;
    }

    for (Index source_id : registry_.AllIds()) {
      if (source_id == target_segment_id) {
        continue;
      }
      if (!registry_.Has(source_id, source_state)) {
        continue;
      }
      const PolarSegment& source_segment =
          registry_.Get(source_id, source_state);
      // Same baseline-wrap fix as AddFieldAt's own slow path -- see
      // that method's own comment (right above its own equivalent
      // computation) for why this exists. Kept in sync with it
      // deliberately, since this dump exists specifically to let a
      // comparison catch exactly this kind of discrepancy; a dump that
      // didn't carry the same fix as the method it's diagnosing would
      // stop being representative of what AddFieldAt actually does.
      const Eigen::Vector3d raw_offset =
          target.getPos() - source_segment.getPos();
      const Eigen::Vector3d frac = box_.inverse() * raw_offset;
      const Eigen::Vector3d wrapped_frac =
          frac - frac.array().round().matrix();
      const Eigen::Vector3d min_image_offset = box_ * wrapped_frac;
      const Eigen::Vector3d baseline_shift = raw_offset - min_image_offset;
      for (Index idx = shell_start; idx < shell_end; ++idx) {
        const Eigen::Vector3d t = baseline_shift + translations_[idx].t;
        const Eigen::Vector3d shifted_pos = source_segment[0].getPos() + t;
        const double r = (target.getPos() - shifted_pos).norm();
        dump << target_segment_id << "," << source_id << "," << r << ","
             << t.x() << "," << t.y() << "," << t.z() << "\n";
      }
    }

    const double shell_radius =
        translations_[shell_end > shell_start ? shell_end - 1 : shell_start]
            .r;
    // See DumpNeighborListAppend's own earlier version's own comment
    // (now folded in here) for why this uses only the shell_radius >=
    // r_min_ half of AddFieldAt's own convergence check.
    if (shell_radius >= r_min_) {
      converged = true;
      break;
    }
    shell_start = shell_end;
  }
  (void)converged;
  dump.close();
}

void EwaldRealSpaceSum::DumpPerPairFieldAppend(
    Index target_segment_id, PolarSite& target, EwaldChargeState source_state,
    const std::string& filename) const {
  // See this method's own declaration for what this is for and why.
  std::ofstream dump(filename, std::ios::app);
  dump.precision(15);

  const Eigen::Vector3d original_V = target.V();

  Index shell_start = 0;
  double shell_edge = 0.0;
  bool converged = false;

  while (shell_start < Index(translations_.size())) {
    shell_edge = translations_[shell_start].r +
                (shell_edge > translations_[shell_start].r ? 0.0
                                                            : shell_width_);
    Index shell_end = shell_start;
    while (shell_end < Index(translations_.size()) &&
          translations_[shell_end].r <= shell_edge) {
      ++shell_end;
    }

    const Eigen::Vector3d before_shell = target.V();

    for (Index source_id : registry_.AllIds()) {
      if (source_id == target_segment_id) {
        continue;
      }
      if (!registry_.Has(source_id, source_state)) {
        continue;
      }
      const PolarSegment& source_segment =
          registry_.Get(source_id, source_state);

      // Same baseline-wrap fix as AddFieldAt's own slow path -- see
      // that method's own comment for why this exists. Kept in sync
      // with it deliberately, for the same reason given in
      // DumpNeighborListAppend's own matching comment.
      const Eigen::Vector3d raw_offset =
          target.getPos() - source_segment.getPos();
      const Eigen::Vector3d frac = box_.inverse() * raw_offset;
      const Eigen::Vector3d wrapped_frac =
          frac - frac.array().round().matrix();
      const Eigen::Vector3d min_image_offset = box_ * wrapped_frac;
      const Eigen::Vector3d baseline_shift = raw_offset - min_image_offset;

      for (Index idx = shell_start; idx < shell_end; ++idx) {
        const Eigen::Vector3d t = baseline_shift + translations_[idx].t;
        Index site_idx = 0;
        for (const PolarSite& source_site : source_segment) {
          PolarSite shifted = source_site;
          shifted.setPos(source_site.getPos() + t);

          const Eigen::Vector3d before_pair = target.V();
          interactor_.ApplyInducedField<Estatic::V>(shifted, target);
          const Eigen::Vector3d pair_field = target.V() - before_pair;

          dump << target_segment_id << "," << source_id << ","
               << site_idx << "," << t.x() << "," << t.y() << "," << t.z()
               << "," << pair_field.x() << "," << pair_field.y() << ","
               << pair_field.z() << "\n";
          ++site_idx;
        }
      }
    }
    (void)before_shell;

    const double shell_radius =
        translations_[shell_end > shell_start ? shell_end - 1 : shell_start]
            .r;
    if (shell_radius >= r_min_) {
      converged = true;
      break;
    }
    shell_start = shell_end;
  }
  (void)converged;

  target.V() = original_V;
  dump.close();
}

void EwaldRealSpaceSum::DumpPerPairStaticFieldAppend(
    Index target_segment_id, PolarSite& target, EwaldChargeState source_state,
    const std::string& filename) const {
  // See this method's own declaration for what this is for and why --
  // identical structure to DumpPerPairFieldAppend, just isolating
  // ApplyStaticField's own contribution instead of ApplyInducedField's.
  std::ofstream dump(filename, std::ios::app);
  dump.precision(15);

  const Eigen::Vector3d original_V = target.V();

  Index shell_start = 0;
  double shell_edge = 0.0;
  bool converged = false;

  while (shell_start < Index(translations_.size())) {
    shell_edge = translations_[shell_start].r +
                (shell_edge > translations_[shell_start].r ? 0.0
                                                            : shell_width_);
    Index shell_end = shell_start;
    while (shell_end < Index(translations_.size()) &&
          translations_[shell_end].r <= shell_edge) {
      ++shell_end;
    }

    const Eigen::Vector3d before_shell = target.V();

    for (Index source_id : registry_.AllIds()) {
      if (source_id == target_segment_id) {
        continue;
      }
      if (!registry_.Has(source_id, source_state)) {
        continue;
      }
      const PolarSegment& source_segment =
          registry_.Get(source_id, source_state);

      const Eigen::Vector3d raw_offset =
          target.getPos() - source_segment.getPos();
      const Eigen::Vector3d frac = box_.inverse() * raw_offset;
      const Eigen::Vector3d wrapped_frac =
          frac - frac.array().round().matrix();
      const Eigen::Vector3d min_image_offset = box_ * wrapped_frac;
      const Eigen::Vector3d baseline_shift = raw_offset - min_image_offset;

      for (Index idx = shell_start; idx < shell_end; ++idx) {
        const Eigen::Vector3d t = baseline_shift + translations_[idx].t;
        Index site_idx = 0;
        for (const PolarSite& source_site : source_segment) {
          PolarSite shifted = source_site;
          shifted.setPos(source_site.getPos() + t);

          const Eigen::Vector3d before_pair = target.V();
          interactor_.ApplyStaticField<PolarSite, Estatic::V>(shifted,
                                                               target);
          const Eigen::Vector3d pair_field = target.V() - before_pair;

          dump << target_segment_id << "," << source_id << ","
               << site_idx << "," << t.x() << "," << t.y() << "," << t.z()
               << "," << pair_field.x() << "," << pair_field.y() << ","
               << pair_field.z() << "\n";
          ++site_idx;
        }
      }
    }
    (void)before_shell;

    const double shell_radius =
        translations_[shell_end > shell_start ? shell_end - 1 : shell_start]
            .r;
    if (shell_radius >= r_min_) {
      converged = true;
      break;
    }
    shell_start = shell_end;
  }
  (void)converged;

  target.V() = original_V;
  dump.close();
}
void EwaldRealSpaceSum::DumpCachedNeighborListAppend(
    Index target_segment_id, const PolarSite& target,
    EwaldChargeState source_state, const std::string& filename) const {
  // See this method's own declaration for what this is for.
  std::ofstream dump(filename, std::ios::app);
  dump.precision(15);

  const std::pair<const PolarSite*, EwaldChargeState> cache_key(&target,
                                                                 source_state);
  auto cached = neighbor_cache_.find(cache_key);
  if (cached == neighbor_cache_.end()) {
    dump << target_segment_id << ",NOT_CACHED,,,,,\n";
    dump.close();
    return;
  }

  for (const auto& entry : cached->second) {
    const Index source_id = std::get<0>(entry);
    const Index translation_idx = std::get<1>(entry);
    const Eigen::Vector3d& baseline_shift = std::get<2>(entry);
    const Eigen::Vector3d t = baseline_shift + translations_[translation_idx].t;
    dump << target_segment_id << "," << source_id << "," << translation_idx
         << "," << t.x() << "," << t.y() << "," << t.z() << ","
         << t.norm() << "\n";
  }
  dump.close();
}

}  // namespace xtp
}  // namespace votca
