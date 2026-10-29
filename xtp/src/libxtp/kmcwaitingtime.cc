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
#include <limits>
#include <stdexcept>

// Local VOTCA includes
#include "votca/xtp/kmcwaitingtime.h"

namespace votca {
namespace xtp {

double TimeDependentRates::EscapeRate(Index site, double t) const {
  double sum = 0.0;
  for (Index e = graph_.EdgeOffset(site); e < graph_.EdgeOffset(site + 1);
       ++e) {
    sum += Rate(e, t);
  }
  return sum;
}

namespace {
// cumulative sums of per-edge values within each site
std::vector<double> CumulativePerSite(const ToFGraph& graph,
                                      const std::vector<double>& value,
                                      std::vector<double>& total) {
  std::vector<double> cumulative(value.size());
  total.assign(graph.NumberOfSites(), 0.0);
  for (Index i = 0; i < graph.NumberOfSites(); ++i) {
    double sum = 0.0;
    for (Index e = graph.EdgeOffset(i); e < graph.EdgeOffset(i + 1); ++e) {
      sum += value[e];
      cumulative[e] = sum;
    }
    total[i] = sum;
  }
  return cumulative;
}

// edge of site i with cumulative value reaching target
Index FindEdge(const ToFGraph& graph, const std::vector<double>& cumulative,
               Index site, double target) {
  Index begin = graph.EdgeOffset(site);
  Index end = graph.EdgeOffset(site + 1);
  Index e = Index(std::lower_bound(cumulative.begin() + begin,
                                   cumulative.begin() + end, target) -
                  cumulative.begin());
  return std::min(e, end - 1);
}

double Infinity() { return std::numeric_limits<double>::infinity(); }
}  // namespace

// exponential

void ExponentialWaitingTime::Prepare(const TimeDependentRates& rates) {
  if (!rates.isStatic()) {
    throw std::runtime_error(
        "ExponentialWaitingTime needs static rates, use thinning or the "
        "integrated-rate inversion for time-dependent rates");
  }
  graph_ = &rates.Graph();
  std::vector<double> value(graph_->NumberOfEdges());
  for (Index e = 0; e < graph_->NumberOfEdges(); ++e) {
    value[e] = rates.Rate(e, 0.0);
  }
  cumulative_ = CumulativePerSite(*graph_, value, escape_);
}

double ExponentialWaitingTime::NextTime(Index site, double t,
                                        KMCRandom& random) const {
  double escape = escape_[site];
  return (escape > 0.0) ? t + random.Exponential(escape) : Infinity();
}

Index ExponentialWaitingTime::ChooseHop(Index site, double,
                                        KMCRandom& random) const {
  return FindEdge(*graph_, cumulative_, site, random.Uniform() * escape_[site]);
}

// thinning

void ThinningWaitingTime::Prepare(const TimeDependentRates& rates) {
  rates_ = &rates;
  const ToFGraph& graph = rates.Graph();
  std::vector<double> bound(graph.NumberOfEdges());
  for (Index e = 0; e < graph.NumberOfEdges(); ++e) {
    bound[e] = rates.UpperBound(e);
    if (!(bound[e] >= 0.0)) {
      throw std::runtime_error("ThinningWaitingTime: invalid upper bound");
    }
  }
  cumulative_ = CumulativePerSite(graph, bound, bound_escape_);
}

double ThinningWaitingTime::NextTime(Index site, double t,
                                     KMCRandom& random) const {
  double escape = bound_escape_[site];
  return (escape > 0.0) ? t + random.Exponential(escape) : Infinity();
}

Index ThinningWaitingTime::ChooseHop(Index site, double t,
                                     KMCRandom& random) const {
  const ToFGraph& graph = rates_->Graph();
  Index e = FindEdge(graph, cumulative_, site,
                     random.Uniform() * bound_escape_[site]);
  double bound = rates_->UpperBound(e);
  double rate = rates_->Rate(e, t);
  if (rate > bound * (1.0 + 1e-12)) {
    throw std::runtime_error(
        "ThinningWaitingTime: rate above its upper bound for edge " +
        std::to_string(e));
  }
  return (random.Uniform() * bound <= rate) ? e : -1;
}

// integrated-rate inversion: template for the implementation

void IntegratedRateInversion::Prepare(const TimeDependentRates& rates) {
  rates_ = &rates;
  if (settings_.quadrature_points < 2 || settings_.interpolation_points < 2) {
    throw std::runtime_error(
        "IntegratedRateInversion: need at least two quadrature and "
        "interpolation points");
  }
  if (!(rates.Period() > 0.0)) {
    throw std::runtime_error(
        "IntegratedRateInversion needs periodic rates (Period() > 0)");
  }
  // TODO: per site i, Lambda_i(t) = int_0^t EscapeRate(i, s) ds on [0, T]
  // (settings_.quadrature_points), I_T = Lambda_i(T), monotone interpolant
  throw std::logic_error(
      "IntegratedRateInversion: not implemented yet (template, see "
      "kmcwaitingtime.h)");
}

double IntegratedRateInversion::NextTime(Index, double, KMCRandom&) const {
  // TODO: g = -ln(u); split t into periods and phase; m and the remainder
  // from I_T; monotone root search on [0, 2T) (BEP eqs. 11-18)
  throw std::logic_error("IntegratedRateInversion: not implemented yet");
}

Index IntegratedRateInversion::ChooseHop(Index, double, KMCRandom&) const {
  // TODO: edge e with probability Rate(e, t) / EscapeRate(i, t)
  throw std::logic_error("IntegratedRateInversion: not implemented yet");
}

std::unique_ptr<WaitingTimeStrategy> MakeWaitingTimeStrategy(
    const std::string& name) {
  if (name == "exponential") {
    return std::make_unique<ExponentialWaitingTime>();
  } else if (name == "thinning") {
    return std::make_unique<ThinningWaitingTime>();
  } else if (name == "inversion") {
    return std::make_unique<IntegratedRateInversion>(
        IntegratedRateInversion::Settings());
  }
  throw std::runtime_error("unknown waiting-time strategy '" + name +
                           "', use exponential, thinning or inversion");
}

}  // namespace xtp
}  // namespace votca
