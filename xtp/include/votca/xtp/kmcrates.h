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
#ifndef VOTCA_XTP_KMCRATES_H
#define VOTCA_XTP_KMCRATES_H

// Standard includes
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "tofgraph.h"

namespace votca {
namespace xtp {

/**
 * \brief Hop rates of the edges of a ToFGraph as functions of time.
 *
 * Edges are identified by their global index in the graph (sites i own the
 * edges [ToFGraph::EdgeOffset(i), ToFGraph::EdgeOffset(i+1))). Times are the
 * absolute simulation time in s (the same clock as the injection times), rates
 * in 1/s. Implementations must be thread safe for the const methods (runs are
 * done in parallel).
 *
 * Used by the waiting-time strategies of the KMC engine (WaitingTimeStrategy):
 * the strategy decides how waiting times are drawn from these rates.
 *
 * Example for a time-dependent implementation (to be written for oscillating
 * fields): Marcus rates with the field F(t) = F0 + F1 cos(2 pi f t + phi),
 * evaluated per edge from constants stored at construction (prefactor,
 * site-energy difference, R_ij, reorganization energy), with Period() = 1/f
 * and UpperBound() = the maximum over one period.
 */
class TimeDependentRates {
 public:
  virtual ~TimeDependentRates() = default;

  /// rate of edge e at time t
  virtual double Rate(Index edge, double t) const = 0;

  /// sum of the rates of the out-edges of site i at time t; override if it
  /// can be computed faster than the sum over the edges
  virtual double EscapeRate(Index site, double t) const;

  /// an upper bound of Rate(edge, t) for all t (needed by thinning)
  virtual double UpperBound(Index edge) const = 0;

  /// common period T of all rates if they are periodic, 0 otherwise
  virtual double Period() const { return 0.0; }

  /// true if no rate depends on time
  virtual bool isStatic() const { return false; }

  const ToFGraph& Graph() const { return graph_; }

 protected:
  explicit TimeDependentRates(const ToFGraph& graph) : graph_(graph) {}
  const ToFGraph& graph_;
};

/// the rates stored in the ToFGraph, constant in time
class StaticRates final : public TimeDependentRates {
 public:
  explicit StaticRates(const ToFGraph& graph) : TimeDependentRates(graph) {}
  double Rate(Index edge, double) const override {
    return graph_.getEdge(edge).rate;
  }
  double EscapeRate(Index site, double) const override {
    return graph_.EscapeRate(site);
  }
  double UpperBound(Index edge) const override {
    return graph_.getEdge(edge).rate;
  }
  bool isStatic() const override { return true; }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_KMCRATES_H
