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
#ifndef VOTCA_XTP_GRAPHRANDOMWALK_H
#define VOTCA_XTP_GRAPHRANDOMWALK_H

// Standard includes
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

// Third party includes
#include <Eigen/SparseCore>

// Local VOTCA includes
#include "eigen.h"
#include "tofgraph.h"

namespace votca {
namespace xtp {

/**
 * \brief Expected time of flight from hitting times of a continuous-time
 * random walk on a ToFGraph, without sampling (graph random walk, Chen, van
 * der Hoorn, Baumeier, arXiv:2405.15836).
 *
 * One carrier: the mean hitting times tau of the sink sites solve
 *   Q tau = 1,   Q = diag(D) - W
 * on the transient sites (all sites except sinks from which a sink can be
 * reached), with D_i the escape rate of site i (hops into sinks included)
 * and W_ij the rate of the hop i -> j. Sinks are not unknowns. The second
 * moment follows from Q m2 = 2 tau, the expected time spent on each site for
 * injection distribution w from the adjoint Q^T g = w.
 *
 * Two carriers with site exclusion: the unknowns are pairs of occupied
 * transient sites i < j, numbered j(j-1)/2 + i. Hops into an occupied site
 * are not allowed (they enter neither D nor W). Two observables:
 *  - First: time until the first carrier reaches a sink (the paper's
 *    definition), Q_PP tau = 1.
 *  - Transit: carriers are removed when they reach a sink, tau is the mean
 *    transit time per carrier, (1/2) E[int N(t) dt]. The pair block is the
 *    same Q_PP; reaching a sink moves to the one-carrier state of the
 *    remaining carrier: Q_PP tau_P = 1 + W_PS tau_S with Q tau_S = 1/2.
 *
 * If the rates obey detailed balance (ToFGraph::isReversible), Q and Q_PP
 * are similar to symmetric positive definite matrices S = P Q P^-1 with
 * P = diag(exp(-U/2kT)) (for pairs the product of the two sites), solved by
 * sparse Cholesky (LDL^T) or, above Settings::direct_limit unknowns, Jacobi
 * preconditioned conjugate gradients; otherwise sparse LU or BiCGSTAB with
 * incomplete LU. Every solution is refined against Q itself.
 *
 * The iterative tolerance has to be tight: with energetic disorder the
 * hitting times are dominated by slow (trap) modes that a relative residual
 * of 1e-10 can leave unresolved (7 % error in tau for a 40^3 lattice with
 * sigma = 0.19 eV), hence the default of 1e-14. Direct factorization of
 * graphs from 3D morphologies costs O(n^2) (O(n^2.5) in the pair space),
 * which is why the iterative solvers take over above a few thousand
 * unknowns.
 *
 * Injection is uniform over the connected source sites (one carrier) or
 * over pairs of them (two carriers). Times are in s.
 */
class GraphRandomWalk {
 public:
  enum class Solver { Auto, Cholesky, LU, CG, BiCGSTAB };
  enum class Observable { First, Transit };

  struct Settings {
    /// auto: direct solvers up to this many unknowns, iterative above
    Index direct_limit = 5000;
    /// relative residual of the iterative solvers
    double tolerance = 1e-14;
    Index max_iterations = 100000;
    /// two carriers: refuse larger pair spaces (memory)
    Index max_states = 50000000;
  };

  struct Result {
    std::string solver;
    Index carriers = 1;
    Observable observable = Observable::First;
    Index unknowns = 0;
    /// largest componentwise backward error |Q x - b|_i /
    /// (|Q| |x| + |b|)_i after refinement, over all solves
    double backward_error = 0.0;
    Index refinement_steps = 0;
    /// iterative solvers: total iterations
    Index iterations = 0;
    /// false if an iterative solve hit the iteration limit or the backward
    /// error after refinement is above 1e-10
    bool converged = true;

    /// injection states: one carrier: one site, two carriers: two sites
    /// (graph indices); the mean time per state and (observable First only)
    /// the standard deviation of the hitting time
    std::vector<std::vector<Index>> sources;
    Eigen::VectorXd tau;
    Eigen::VectorXd tau_std;

    /// uniform injection: mean, and (First only) standard deviation of the
    /// hitting time
    double tau_mean = 0.0;
    double tau_std_injection = 0.0;
    /// harmonic mean over the injection states (eq. 3.4 of the paper)
    double tau_harmonic = 0.0;

    /// expected time each graph site is occupied before absorption, uniform
    /// injection; sums to carriers x tau_mean for Transit and for one
    /// carrier, to 2 x tau_mean for two carriers and First. Zero on sinks.
    Eigen::VectorXd occupation_time;
  };

  explicit GraphRandomWalk(const ToFGraph& graph) : graph_(graph) {
    if (graph_.isBulk()) {
      throw std::runtime_error(
          "GraphRandomWalk: hitting times need a ToF graph with sinks, not a "
          "bulk graph without sinks");
    }
  }

  void setSettings(const Settings& settings) { settings_ = settings; }

  Result SingleCarrier(Solver solver = Solver::Auto) const;

  Result TwoCarriers(Observable observable, Solver solver = Solver::Auto) const;

  static Solver SolverFromString(const std::string& name);
  static Observable ObservableFromString(const std::string& name);

 private:
  struct Transient {
    std::vector<Index> index_of;      // graph site -> transient index or -1
    std::vector<Index> site_of;       // transient index -> graph site
    Eigen::VectorXd half_log_weight;  // -(U - Uref)/2kT
  };
  Transient FindTransient() const;

  const ToFGraph& graph_;
  Settings settings_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GRAPHRANDOMWALK_H
