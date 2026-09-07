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
    for (const auto& pair : cached->second) {
      const Index source_id = pair.first;
      const Index translation_idx = pair.second;
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
      const Eigen::Vector3d& t = translations_[translation_idx].t;
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
  std::vector<std::pair<Index, Index>> visited_pairs;

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

      for (Index idx = shell_start; idx < shell_end; ++idx) {
        visited_pairs.emplace_back(source_id, idx);
        const Eigen::Vector3d& t = translations_[idx].t;
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

}  // namespace xtp
}  // namespace votca
