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
#ifndef VOTCA_XTP_TOFSIMULATION_H
#define VOTCA_XTP_TOFSIMULATION_H

// Standard includes
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// Local VOTCA includes
#include "kmcengine.h"
#include "kmcrates.h"
#include "kmcwaitingtime.h"
#include "tofgraph.h"

namespace votca {
namespace xtp {

/**
 * \brief Many independent time-of-flight KMC runs on a ToFGraph and their
 * statistics.
 *
 * Per run: the time until the first carrier reaches a sink ("first") and
 * the mean transit time per absorbed carrier, transit = absorption -
 * injection time ("transit"); with all carriers injected at t = 0 these are
 * the observables of GraphRandomWalk. Censored runs (stopped by a limit)
 * are counted and left out of the means, which biases them low; the caller
 * should report them.
 *
 * Transient current (Ramo-Shockley): every accepted hop contributes its
 * displacement along the transport normal, binned on a logarithmic time
 * axis and normalized by L x injected carriers x bin width, so that its
 * integral is the mean displacement per injected carrier in units of L
 * (about 1 for complete collection). The collected fraction is the number
 * of absorptions up to the end of each bin per injected carrier.
 *
 * Runs are done in blocks in parallel (OpenMP); per-run data are reduced
 * in run order, so results do not depend on the number of threads.
 */
class ToFSimulation {
 public:
  struct Settings {
    std::uint64_t seed = 1;
    Index runs = 1000;
    std::vector<KMCEngine::Pulse> pulses = {KMCEngine::Pulse{}};
    KMCEngine::Limits limits;
    double tmin = 1e-15;  // s, current/collection time axis
    double tmax = 1.0;
    Index bins_per_decade = 10;
  };

  struct RunData {
    KMCEngine::RunStatus status;
    double first_arrival = 0.0;  // valid if status.absorbed > 0
    double mean_transit = 0.0;   // valid if status.absorbed > 0
  };

  struct Statistic {
    Index samples = 0;
    double mean = 0.0;
    double standard_error = 0.0;
    double standard_deviation = 0.0;
    /// share of the sum contributed by the largest 1 % of the samples
    double top_percent_share = 0.0;
  };

  struct Result {
    std::vector<RunData> runs;
    Index censored = 0;
    Index injected = 0;
    Index blocked_injections = 0;
    Index hops = 0;
    Index null_events = 0;
    Statistic first;                // over uncensored runs
    Statistic transit;              // over uncensored runs
    std::vector<double> bin_edges;  // s, size bins + 1
    std::vector<double> current;    // 1/s, per bin
    std::vector<double> collected;  // fraction absorbed by the bin end
    /// displacement (in units of L per injected carrier) before tmin and
    /// after tmax
    double displacement_before = 0.0;
    double displacement_after = 0.0;
  };

  /// static rates of the graph, exponential waiting times
  ToFSimulation(const ToFGraph& graph, const Settings& settings);

  /// any rates and waiting-time strategy; calls strategy.Prepare(rates);
  /// rates and strategy must outlive the simulation
  ToFSimulation(const ToFGraph& graph, const Settings& settings,
                const TimeDependentRates& rates, WaitingTimeStrategy& strategy);

  const std::string& StrategyName() const { return strategy_name_; }

  Result Simulate() const;

  static Statistic Analyze(std::vector<double> samples);

 private:
  RunData SingleRun(Index run, std::vector<double>& current,
                    std::vector<double>& arrivals, double& before,
                    double& after) const;
  Index Bin(double t) const;

  void Setup();

  const ToFGraph& graph_;
  Settings settings_;
  std::unique_ptr<StaticRates> own_rates_;
  std::unique_ptr<WaitingTimeStrategy> own_strategy_;
  std::unique_ptr<KMCEngine> engine_;
  std::string strategy_name_;
  std::vector<double> edges_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TOFSIMULATION_H
