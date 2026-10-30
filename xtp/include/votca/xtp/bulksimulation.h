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
#ifndef VOTCA_XTP_BULKSIMULATION_H
#define VOTCA_XTP_BULKSIMULATION_H

// Standard includes
#include <cstdint>
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
 * \brief Kinetic Monte Carlo of a fixed number of carriers with site
 * exclusion in a periodic box (bulk transport, calculator kmcmultiple).
 *
 * Each run places the carriers on distinct random injection sites at t = 0
 * and propagates them with KMCEngine on a bulk ToFGraph (ToFGraph::BulkOptions)
 * until a limit (simulated time, events or wall time). Observables are taken
 * over the window [equilibration, end of run]:
 * - drift velocity: unwrapped displacement per carrier over the window;
 * - diffusion tensor: the window is cut into J blocks of equal length tau
 *   (blocks <= J < 2 blocks); the increments d of all carriers over the
 *   blocks give D = Cov(d) / (2 tau), with the mean (drift) removed. This
 *   is the long-time D if tau is longer than the memory of the motion
 *   (trapping times); compare runs with different numbers of blocks.
 *   Block boundaries are found during the run by doubling tau, so the
 *   positions at the boundaries are exact also when the length of the run
 *   is set by the number of events;
 * - occupation probability of every site (summed over carriers, so the sum
 *   over sites is the number of carriers).
 *
 * Several runs are independent (random numbers from (seed, run)), done in
 * parallel and reduced in run order; with more than one run the standard
 * errors of velocity and D are taken over the runs, otherwise the
 * velocity error comes from the block increments.
 */
class BulkSimulation {
 public:
  struct Settings {
    std::uint64_t seed = 1;
    Index runs = 1;
    Index carriers = 1;
    KMCEngine::Limits limits;
    double equilibration = 0.0;  // s
    Index blocks = 50;
    /// trajectory of run 0 sampled every interval (s) and/or after every
    /// sample_hops accepted hops; 0: none
    double sample_interval = 0.0;
    Index sample_hops = 0;
  };

  /// state of run 0 at a sampling time
  struct Sample {
    double time = 0.0;
    Index hops = 0;
    /// unwrapped positions of the carriers
    std::vector<Eigen::Vector3d> positions;
    /// mean site energy of the carriers
    double energy = 0.0;
  };

  struct RunData {
    KMCEngine::RunStatus status;
    /// observation window [start, end]
    double start = 0.0;
    double end = 0.0;
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    Eigen::Matrix3d diffusion = Eigen::Matrix3d::Zero();
    /// covariance of the block increments of the velocity estimate
    Eigen::Matrix3d velocity_covariance = Eigen::Matrix3d::Zero();
    Index blocks = 0;
    double block_length = 0.0;
    /// displacement of each carrier over the window
    std::vector<Eigen::Vector3d> displacement;
  };

  struct Result {
    std::vector<RunData> runs;
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    Eigen::Vector3d velocity_error = Eigen::Vector3d::Zero();
    Eigen::Matrix3d diffusion = Eigen::Matrix3d::Zero();
    /// standard error over runs, zero for one run
    Eigen::Matrix3d diffusion_error = Eigen::Matrix3d::Zero();
    /// occupation probability per site, averaged over runs
    Eigen::VectorXd occupation;
    std::vector<Sample> trajectory;
    Index hops = 0;
    Index null_events = 0;
    Index wall_time_stops = 0;
  };

  /// one carrier without sampling (master equation): stationary
  /// distribution pi (L^T pi = 0), drift v = sum_i pi_i b_i with
  /// b_i = sum_j w_ij R_ij, and the diffusion tensor from the corrector chi
  /// (L chi = v - b): D = 1/2 sum_i pi_i sum_j w_ij g g^T, g = R_ij + chi_j
  /// - chi_i. Needs a unique stationary distribution: exactly one closed
  /// class of sites (sites from which the carrier cannot leave); otherwise
  /// available is false.
  struct ExactResult {
    bool available = false;
    std::string message;
    Index sites = 0;  // in the closed class
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    Eigen::Matrix3d diffusion = Eigen::Matrix3d::Zero();
    Eigen::VectorXd occupation;  // pi, zero outside the closed class
    std::string solver;
    double backward_error = 0.0;
    bool converged = false;
  };
  /// direct sparse LU up to direct_limit sites, BiCGSTAB with ILUT above
  static ExactResult Exact(const ToFGraph& graph, Index direct_limit = 20000,
                           double tolerance = 1e-12);

  /// graph must be a bulk graph; strategy is prepared for the static rates of
  /// the graph unless a prepared strategy is given
  BulkSimulation(const ToFGraph& graph, const Settings& settings);
  BulkSimulation(const ToFGraph& graph, const Settings& settings,
                 const TimeDependentRates& rates,
                 WaitingTimeStrategy& strategy);

  Result Simulate() const;

 private:
  const ToFGraph& graph_;
  Settings settings_;
  std::unique_ptr<TimeDependentRates> own_rates_;
  std::unique_ptr<WaitingTimeStrategy> own_strategy_;
  std::unique_ptr<KMCEngine> engine_;

  void Check() const;
  RunData SingleRun(Index run, Eigen::VectorXd& occupation,
                    std::vector<Sample>* trajectory) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BULKSIMULATION_H
