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
#ifndef VOTCA_XTP_TRAPFINDER_H
#define VOTCA_XTP_TRAPFINDER_H

// Standard includes
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "tofgraph.h"

namespace votca {
namespace xtp {

/**
 * \brief Identifies traps (single sites or clusters) of a charge transport
 * network from the slow modes of the master equation, after Chen, van der
 * Hoorn and Baumeier (spectral clustering of molecular charge transport
 * networks).
 *
 * The generator -Q = D - W (D: escape rates) of the hopping graph is
 * diagonalized (dense, all sites with hops). Its eigenvalues 0 = l_0 <= l_1
 * <= ... are relaxation rates. A trap shows up as a slow mode separated by
 * a gap from the faster ones: the number of slow modes is the index i
 * (1..modes-2) with the largest ratio l_{i+1}/l_i. The density modes
 * (left eigenvectors) of
 * the slow modes are localized on the traps: the sites with
 * |v_k|/max|v_k| > `threshold` form the trap candidates of mode k; their
 * connected components (in the hopping graph) are traps if they are
 * localized (at most `max_size` sites carrying at least the fraction `mass`
 * of the weight of their side of the mode), and components contained in
 * others are dropped. A trap must hold the carrier: its release rate (see
 * below) times `gap` must not exceed the first fast relaxation rate
 * l_{i+1}. The eigenvalue ratio alone does not separate traps from
 * Gaussian disorder in small networks: when a trap holds most of the
 * equilibrium weight, the slow mode is its capture, which is fast in a
 * small box; the release rate is the discriminating quantity. Unlike the paper,
 * which uses the eigenvectors of the jump chain I - D^-1 W and a subsequent
 * K-means step, the supports of the density modes of the continuous-time
 * generator are used directly: they are localized on the traps also when the
 * trap holds most of the equilibrium weight. A density mode moves weight
 * between the trap and the rest of the network (it sums to zero); the trap
 * is taken on the side with fewer sites.
 *
 * For every trap the equilibrium weight pi(C) (stationary distribution of
 * the network, also without detailed balance), the relaxation rate of its
 * mode, split into release and capture rates in a two-state picture, and
 * the hop rate
 * sum_{i in C, j not in C} pi_i w_ij / pi(C) are given.
 */
class TrapFinder {
 public:
  struct Settings {
    double threshold = 0.9;  // the paper's a
    /// a trap releases a carrier at least gap times slower than the first
    /// fast mode relaxes
    double gap = 10.0;
    /// a trap carries at least this fraction of the weight of its side of
    /// the density mode, on at most max_size sites
    double mass = 0.5;
    Index max_size = 20;
    Index modes = 50;        // eigenvalues considered
    Index max_sites = 5000;  // dense eigenvalue problem
  };

  struct Trap {
    std::vector<Index> sites;  // graph site indices, deepest first
    double weight = 0.0;       // equilibrium probability
    /// relaxation rate of the slow mode that found the trap (1/s); in the
    /// two-state picture of the mode it is release + capture rate, and the
    /// weight is capture / (release + capture):
    double relaxation_rate = 0.0;
    double release_rate = 0.0;  // relaxation_rate (1 - weight)
    double capture_rate = 0.0;  // relaxation_rate weight
    /// equilibrium flux of hops out of the trap per occupation (1/s); much
    /// larger than release_rate when the carrier mostly hops back
    double hop_rate = 0.0;
  };

  struct Result {
    /// the smallest eigenvalues of -Q (relaxation rates, 1/s)
    std::vector<double> eigenvalues;
    double gap = 0.0;                // largest ratio l_{i+1}/l_i
    Index slow_modes = 0;            // i of the largest ratio
    std::vector<Trap> traps;         // empty if gap < Settings::gap
    std::vector<double> occupation;  // stationary distribution per site
  };

  static Result Find(const ToFGraph& graph, const Settings& settings);
  static Result Find(const ToFGraph& graph) { return Find(graph, Settings()); }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TRAPFINDER_H
