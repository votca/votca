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
#ifndef VOTCA_XTP_KMCENGINE_H
#define VOTCA_XTP_KMCENGINE_H

// Standard includes
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "kmcrandom.h"
#include "kmcwaitingtime.h"
#include "tofgraph.h"

namespace votca {
namespace xtp {

/**
 * \brief Kinetic Monte Carlo of carriers with site exclusion on a ToFGraph.
 *
 * Every carrier has its own clock: the pluggable WaitingTimeStrategy gives
 * the time of its next candidate event and, at that time, the hop
 * (KMCWaitingTime.h has the exactness argument). The engine processes the
 * earliest pending event (a min-heap over carriers, merged with the
 * scheduled injections); a candidate hop into an occupied site, or a
 * strategy returning -1, is a null event: time has advanced, the carrier
 * draws a new candidate time. Exclusion therefore never shortens waiting
 * times. Only the clock of the acting carrier changes at an event.
 *
 * Strategies: ExponentialWaitingTime (static rates), ThinningWaitingTime
 * (time-dependent rates with upper bounds), and the template
 * IntegratedRateInversion (periodic rates).
 *
 * Carriers enter at scheduled injection times on a free connected source
 * site (uniform) and leave at a sink. Optionally they also decay with a
 * per-site rate (setDecay): a second, exponential clock competing with the
 * hops, redrawn after every event of the carrier (exact, as the decay is
 * memoryless); a decayed carrier can be replaced at once by a new one on a
 * free source site (constant number of carriers). Observers receive every
 * injection, accepted hop, absorption and decay.
 *
 * Runs are independent and reproducible: the random numbers of a run depend
 * only on (seed, run index), not on threads or the order of runs.
 */
class KMCEngine {
 public:
  using Random = KMCRandom;

  struct Pulse {
    enum class Shape { Delta, Gaussian, Rectangular };
    double start = 0.0;  // s
    Index carriers = 1;
    Shape shape = Shape::Delta;
    double width = 0.0;  // s: standard deviation (Gaussian) or full width
  };

  /// injection times of all carriers of the pulses, sorted; finite-width
  /// pulses draw the times of their carriers from the shape (Gaussian
  /// truncated to t >= 0)
  static std::vector<double> InjectionTimes(const std::vector<Pulse>& pulses,
                                            Random& random);
  static Pulse::Shape ShapeFromString(const std::string& name);

  class Observer {
   public:
    virtual ~Observer() = default;
    virtual void Inject(Index /*carrier*/, Index /*site*/, double /*t*/) {}
    virtual void Hop(Index /*carrier*/, Index /*from*/, Index /*to*/,
                     const Eigen::Vector3d& /*dR*/, double /*t*/) {}
    virtual void Absorb(Index /*carrier*/, Index /*site*/, double /*t*/) {}
    virtual void Decay(Index /*carrier*/, Index /*site*/, double /*t*/) {}
  };

  struct DecaySettings {
    /// decay rate per site (1/s), indexed like the graph sites; nullptr: no
    /// decay. Must outlive the engine.
    const std::vector<double>* rates = nullptr;
    /// inject a new carrier at the time of every decay
    bool replace = false;
  };

  struct Limits {
    /// simulated time, s
    double max_time = std::numeric_limits<double>::infinity();
    /// candidate events (hops and null events)
    Index max_events = std::numeric_limits<Index>::max();
    /// wall-clock time of the run, s (checked every 4096 events)
    double max_wall_time = std::numeric_limits<double>::infinity();
    /// stop after this many decays
    Index max_decays = std::numeric_limits<Index>::max();
  };

  struct RunStatus {
    double end_time = 0.0;
    Index hops = 0;
    /// candidate events without a hop: destination occupied or rejected by
    /// the strategy (thinning)
    Index null_events = 0;
    Index injected = 0;
    Index absorbed = 0;
    Index decayed = 0;
    /// injections dropped because all source sites were occupied
    Index blocked_injections = 0;
    /// stopped by a limit (or no carrier could move) before all carriers
    /// were absorbed
    bool censored = false;
    /// stopped by Limits::max_time (else by events, wall time or because
    /// no carrier could move)
    bool reached_max_time = false;
    bool reached_wall_time = false;
    bool reached_max_decays = false;
    /// candidate events processed
    Index events = 0;
  };

  /// strategy must be prepared (WaitingTimeStrategy::Prepare) for rates on
  /// this graph and outlive the engine
  KMCEngine(const ToFGraph& graph, const WaitingTimeStrategy& strategy)
      : graph_(graph), strategy_(strategy) {}

  /// per-site decay; rates must be >= 0 and finite, one per site
  void setDecay(const DecaySettings& decay);

  /// one run; injection_times sorted; carriers are numbered in order of
  /// injection (scheduled ones and replacements, including injections that
  /// were blocked)
  RunStatus Run(const std::vector<double>& injection_times, Random& random,
                Observer& observer, const Limits& limits) const;

 private:
  const ToFGraph& graph_;
  const WaitingTimeStrategy& strategy_;
  DecaySettings decay_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_KMCENGINE_H
