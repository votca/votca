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
#ifndef VOTCA_XTP_EWALDREGISTRY_H
#define VOTCA_XTP_EWALDREGISTRY_H

// Standard includes
#include <map>
#include <vector>

// Local VOTCA includes
#include "checkpoint.h"
#include "classicalsegment.h"

/**
 * \brief Central, index-based store of PolarSegment objects, keyed by
 *        segment id and charge state.
 *
 * This is the modern replacement for the legacy Ewald code's PolarTop /
 * PolarSeg / APolarSite machinery. Where the legacy design shared raw
 * PolarSeg* pointers across multiple containers (BGN/MGN/FGN/FGC,
 * QM0/MM1/MM2) and reassigned them at runtime with manually-maintained
 * ownership flags, EwaldRegistry instead holds every PolarSegment by value,
 * in one place, and everything else (background machinery, regions, ...)
 * refers to a segment by its stable (id, charge state) key rather than by
 * pointer. "Moving" a segment between roles becomes a change in which key
 * is queried, never a change in who owns the underlying data.
 *
 * A segment that only ever needs one representation (the common case: the
 * neutral background) is registered once, under EwaldChargeState::Neutral.
 * A segment that additionally needs an excited-state representation (the
 * foreground molecule under study) is registered a second time, under
 * EwaldChargeState::Electron or EwaldChargeState::Hole, alongside its
 * existing neutral entry. Both entries remain independently addressable;
 * neither is preferred or "active" by default. Callers explicitly say which
 * state they want when they query.
 *
 * EwaldRegistry deliberately does not know how to build a PolarSegment from
 * a Topology; that is a separate concern (see SegmentMapper). This class is
 * only ever a store.
 */

namespace votca {
namespace xtp {

enum class EwaldChargeState { Neutral, Electron, Hole };

class EwaldRegistry {
 public:
  EwaldRegistry() = default;

  // Registers `segment` under (id, state). Overwrites any existing entry
  // for the same (id, state) pair. `id` is expected to match the stable
  // segment id used elsewhere in the topology (Segment::getId()); it is
  // not required to match segment.getId() and is kept as an explicit
  // argument so the same PolarSegment content could in principle be
  // registered under a different id if that is ever useful, though in
  // normal use the two will agree.
  void Register(Index id, EwaldChargeState state, PolarSegment segment);

  // True if an entry exists for (id, state).
  bool Has(Index id, EwaldChargeState state) const;

  // Returns the stored segment for (id, state). Throws std::out_of_range
  // if no such entry exists; callers are expected to check Has() first
  // when the presence of a given state is not already guaranteed by
  // construction.
  const PolarSegment& Get(Index id, EwaldChargeState state) const;
  PolarSegment& Get(Index id, EwaldChargeState state);

  // Removes the entry for (id, state), if present. No-op otherwise.
  void Erase(Index id, EwaldChargeState state);

  // Every segment id that has at least one registered state.
  std::vector<Index> AllIds() const;

  // Every state currently registered for `id`. Empty if `id` is unknown.
  std::vector<EwaldChargeState> StatesFor(Index id) const;

  Index size() const { return static_cast<Index>(store_.size()); }

  void WriteToCpt(CheckpointWriter& w) const;
  void ReadFromCpt(CheckpointReader& r);

 private:
  using Key = std::pair<Index, EwaldChargeState>;

  struct KeyLess {
    bool operator()(const Key& a, const Key& b) const {
      if (a.first != b.first) {
        return a.first < b.first;
      }
      return static_cast<int>(a.second) < static_cast<int>(b.second);
    }
  };

  std::map<Key, PolarSegment, KeyLess> store_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDREGISTRY_H
