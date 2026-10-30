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
#ifndef VOTCA_XTP_LIFETIMESIMULATION_H
#define VOTCA_XTP_LIFETIMESIMULATION_H

// Standard includes
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "kmcengine.h"
#include "kmcrates.h"
#include "kmcwaitingtime.h"
#include "tofgraph.h"

namespace votca {
namespace xtp {

/**
 * \brief Kinetic Monte Carlo of carriers (excitons) that hop and decay with
 * a site-dependent rate, at constant number: every decayed carrier is
 * replaced at once by a new one on a random free injection site (calculator
 * kmclifetime).
 *
 * Each run starts with the carriers on distinct random injection sites and
 * stops after its share of the requested decays. Per decay: lifetime
 * (decay - injection time), displacement (unwrapped), hops and decay site.
 * Statistics are over all decays; their standard errors assume independent
 * decays (exact for one carrier). Also the occupation probability of every
 * site (summed over carriers) and the fraction of decays per site. Runs are
 * done in parallel and reduced in run order.
 *
 * Exact(): for one carrier the same averages without sampling, from linear
 * systems with A = K + diag(W) - W (K: decay rates, W: hop rates):
 * A tau = 1 (lifetime), A m = b with b_i = sum_j w_ij R_ij (mean
 * displacement), A s = sum_j w_ij (R_ij^2 + 2 R_ij m_j) (mean squared
 * displacement, per component), and A^T g = p0 for the expected time per
 * site (occupation) and the decay probabilities k_j g_j, all averaged over a
 * uniform start on the injection sites.
 */
class LifetimeSimulation {
 public:
  struct Settings {
    std::uint64_t seed = 1;
    Index runs = 1;
    Index carriers = 1;
    /// decays in total, distributed over the runs
    Index decays = 1000;
    /// wall-clock limit of each run (s)
    double max_wall_time = std::numeric_limits<double>::infinity();
    /// keep the record of every decay
    bool keep_records = false;
    /// energy trace of run 0: mean site energy of the carriers after every
    /// energy_hops hops (0: none)
    Index energy_hops = 0;
  };

  struct Record {
    Index run = 0;
    Index carrier = 0;
    double time = 0.0;      // of the decay
    double lifetime = 0.0;  // decay - injection
    Index hops = 0;
    Index site = 0;
    Eigen::Vector3d displacement = Eigen::Vector3d::Zero();
  };

  /// sums over the decays of one run
  struct DecaySums {
    Index decays = 0;
    double lifetime = 0.0;
    double lifetime2 = 0.0;
    double distance = 0.0;
    double distance2 = 0.0;  // |R|^2
    double distance4 = 0.0;  // |R|^4
    Eigen::Vector3d components = Eigen::Vector3d::Zero();
    double hops = 0.0;
  };

  struct EnergySample {
    double time = 0.0;
    Index hops = 0;
    double energy = 0.0;  // mean site energy of the carriers
  };

  struct Result {
    Index decays = 0;
    double lifetime = 0.0;
    double lifetime_error = 0.0;
    double hops = 0.0;           // mean per carrier
    double mean_distance = 0.0;  // <|R|>
    double mean_distance_error = 0.0;
    double squared_displacement = 0.0;  // <|R|^2>
    double squared_displacement_error = 0.0;
    Eigen::Vector3d squared_components = Eigen::Vector3d::Zero();  // <R_a^2>
    /// fraction of the decays per site
    Eigen::VectorXd decay_fraction;
    /// occupation probability per site (sum: number of carriers)
    Eigen::VectorXd occupation;
    std::vector<Record> records;
    std::vector<EnergySample> energy;
    Index hop_events = 0;
    Index null_events = 0;
    Index wall_time_stops = 0;
    double simulated_time = 0.0;  // all runs
  };

  struct ExactResult {
    double lifetime = 0.0;
    double squared_displacement = 0.0;
    Eigen::Vector3d squared_components = Eigen::Vector3d::Zero();
    Eigen::VectorXd decay_probability;
    /// expected time per site per carrier, normalized by the lifetime
    Eigen::VectorXd occupation;
    std::string solver;
    double backward_error = 0.0;
    bool converged = false;
  };

  /// graph: a bulk graph (ToFGraph::BulkOptions) for the usual setup; decay
  /// rates per site (1/s), positive for every site that takes part
  LifetimeSimulation(const ToFGraph& graph, std::vector<double> decay_rates,
                     const Settings& settings);
  LifetimeSimulation(const ToFGraph& graph, std::vector<double> decay_rates,
                     const Settings& settings, const TimeDependentRates& rates,
                     WaitingTimeStrategy& strategy);

  Result Simulate() const;

  /// one carrier, static rates of the graph; direct sparse LU up to
  /// direct_limit sites, BiCGSTAB with ILUT above
  static ExactResult Exact(const ToFGraph& graph,
                           const std::vector<double>& decay_rates,
                           Index direct_limit = 20000,
                           double tolerance = 1e-12);

 private:
  const ToFGraph& graph_;
  std::vector<double> decay_rates_;
  Settings settings_;
  std::unique_ptr<TimeDependentRates> own_rates_;
  std::unique_ptr<WaitingTimeStrategy> own_strategy_;
  std::unique_ptr<KMCEngine> engine_;

  void Check() const;
  void Setup(WaitingTimeStrategy& strategy, const TimeDependentRates& rates);
  struct RunData {
    KMCEngine::RunStatus status;
    std::vector<Record> records;
    Eigen::VectorXd occupation;  // time per site
    Eigen::VectorXd decays;      // per site
    DecaySums sums;
  };
  RunData SingleRun(Index run, Index decays,
                    std::vector<EnergySample>* energy) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_LIFETIMESIMULATION_H
