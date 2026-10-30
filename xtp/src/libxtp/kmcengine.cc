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
#include <chrono>
#include <cmath>
#include <functional>
#include <iterator>
#include <queue>
#include <stdexcept>
#include <utility>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/kmcengine.h"

namespace votca {
namespace xtp {

KMCRandom::KMCRandom(std::uint64_t seed, std::uint64_t stream) {
  // std::seed_seq and mt19937_64 are fully specified by the standard
  std::seed_seq seq{
      std::uint32_t(seed & 0xffffffffu), std::uint32_t(seed >> 32),
      std::uint32_t(stream & 0xffffffffu), std::uint32_t(stream >> 32)};
  engine_.seed(seq);
}

Index KMCRandom::Integer(Index n) {
  if (n <= 0) {
    throw std::runtime_error("KMCRandom::Integer: n must be positive");
  }
  const std::uint64_t un = std::uint64_t(n);
  const std::uint64_t limit = std::numeric_limits<std::uint64_t>::max() -
                              std::numeric_limits<std::uint64_t>::max() % un;
  std::uint64_t x;
  do {
    x = engine_();
  } while (x >= limit);
  return Index(x % un);
}

double KMCRandom::Normal() {
  double u1 = Uniform();
  double u2 = Uniform();
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * tools::conv::Pi * u2);
}

KMCEngine::Pulse::Shape KMCEngine::ShapeFromString(const std::string& name) {
  if (name == "delta") {
    return Pulse::Shape::Delta;
  } else if (name == "gaussian") {
    return Pulse::Shape::Gaussian;
  } else if (name == "rectangular") {
    return Pulse::Shape::Rectangular;
  }
  throw std::runtime_error("KMCEngine: unknown pulse shape '" + name +
                           "', use delta, gaussian or rectangular");
}

std::vector<double> KMCEngine::InjectionTimes(const std::vector<Pulse>& pulses,
                                              Random& random) {
  std::vector<double> times;
  for (const Pulse& p : pulses) {
    if (p.carriers < 0 || p.start < 0.0 || p.width < 0.0) {
      throw std::runtime_error(
          "KMCEngine: pulses need start >= 0, width >= 0 and carriers >= 0");
    }
    for (Index c = 0; c < p.carriers; ++c) {
      double t = p.start;
      if (p.shape == Pulse::Shape::Gaussian && p.width > 0.0) {
        do {
          t = p.start + p.width * random.Normal();
        } while (t < 0.0);
      } else if (p.shape == Pulse::Shape::Rectangular) {
        t = p.start + p.width * random.Uniform();
      }
      times.push_back(t);
    }
  }
  std::sort(times.begin(), times.end());
  return times;
}

void KMCEngine::setDecay(const DecaySettings& decay) {
  if (decay.rates != nullptr) {
    if (Index(decay.rates->size()) != graph_.NumberOfSites()) {
      throw std::runtime_error(
          "KMCEngine: need one decay rate per site of the graph");
    }
    if (!std::all_of(decay.rates->begin(), decay.rates->end(),
                     [](double k) { return k >= 0.0 && std::isfinite(k); })) {
      throw std::runtime_error(
          "KMCEngine: decay rates must be finite and not negative");
    }
  } else if (decay.replace) {
    throw std::runtime_error("KMCEngine: replacement needs decay rates");
  }
  decay_ = decay;
}

KMCEngine::RunStatus KMCEngine::Run(const std::vector<double>& injection_times,
                                    Random& random, Observer& observer,
                                    const Limits& limits) const {
  if (!std::is_sorted(injection_times.begin(), injection_times.end())) {
    throw std::runtime_error("KMCEngine: injection times must be sorted");
  }
  const double inf = std::numeric_limits<double>::infinity();
  const Index n = graph_.NumberOfSites();
  const Index ncarriers = Index(injection_times.size());
  const std::vector<Index>& sources = graph_.SourceSites();
  std::vector<Index> occupant(n, -1);
  // per carrier, indexed by carrier number (replacements add carriers)
  std::vector<Index> site_of(ncarriers, -1);
  std::vector<char> decays_next(ncarriers, 0);
  const std::vector<double>* decay = decay_.rates;
  // pending candidate event per carrier: (time, carrier), earliest first;
  // ties broken by carrier index for reproducibility
  using Pending = std::pair<double, Index>;
  std::priority_queue<Pending, std::vector<Pending>, std::greater<Pending>>
      pending;
  RunStatus status;
  Index next = 0;  // next scheduled injection
  Index next_carrier = 0;
  Index alive = 0;
  double t = 0.0;
  Index events = 0;
  const auto wall_start = std::chrono::steady_clock::now();
  auto wall_exceeded = [&]() {
    if ((events & 4095) != 0 || limits.max_wall_time == inf) {
      return false;
    }
    std::chrono::duration<double> elapsed =
        std::chrono::steady_clock::now() - wall_start;
    return elapsed.count() > limits.max_wall_time;
  };

  auto schedule = [&](Index carrier) {
    double tc = strategy_.NextTime(site_of[carrier], t, random);
    // t + dt can round to t when dt < t * 1e-16 (rare short waiting times in
    // long runs); an event at the same time is fine, an earlier one is not
    if (!(tc >= t)) {
      throw std::runtime_error(
          "KMCEngine: the waiting-time strategy returned a time before the "
          "current time");
    }
    decays_next[carrier] = 0;
    if (decay != nullptr) {
      // competing decay clock, memoryless: redrawn after every event
      double k = (*decay)[site_of[carrier]];
      if (k > 0.0) {
        double td = t + random.Exponential(k);
        if (td < tc) {
          tc = td;
          decays_next[carrier] = 1;
        }
      }
    }
    if (tc < inf) {
      pending.emplace(tc, carrier);
    }
  };

  auto inject = [&]() {
    Index carrier = next_carrier++;
    if (carrier >= Index(site_of.size())) {
      site_of.resize(carrier + 1, -1);
      decays_next.resize(carrier + 1, 0);
    }
    std::vector<Index> free;
    free.reserve(sources.size());
    std::copy_if(sources.begin(), sources.end(), std::back_inserter(free),
                 [&occupant](Index s) { return occupant[s] < 0; });
    if (free.empty()) {
      ++status.blocked_injections;
      return;
    }
    Index s = free[random.Integer(Index(free.size()))];
    occupant[s] = carrier;
    site_of[carrier] = s;
    ++alive;
    ++status.injected;
    observer.Inject(carrier, s, t);
    schedule(carrier);
  };

  while (alive > 0 || next < ncarriers) {
    double t_inject = (next < ncarriers) ? injection_times[next] : inf;
    double t_event = pending.empty() ? inf : pending.top().first;
    if (t_inject <= t_event) {
      if (t_inject == inf || t_inject > limits.max_time) {
        status.censored = true;  // nobody can move any more, or time limit
        status.reached_max_time = t_inject != inf;
        break;
      }
      t = t_inject;
      ++next;
      inject();
      continue;
    }
    if (t_event > limits.max_time || events >= limits.max_events) {
      status.censored = true;
      status.reached_max_time = t_event > limits.max_time;
      break;
    }
    if (wall_exceeded()) {
      status.censored = true;
      status.reached_wall_time = true;
      break;
    }
    Index carrier = pending.top().second;
    pending.pop();
    t = t_event;
    ++events;
    Index s = site_of[carrier];
    if (decays_next[carrier]) {
      occupant[s] = -1;
      site_of[carrier] = -1;
      --alive;
      ++status.decayed;
      observer.Decay(carrier, s, t);
      if (status.decayed >= limits.max_decays) {
        status.reached_max_decays = true;
        break;
      }
      if (decay_.replace) {
        inject();
      }
      continue;
    }
    Index e = strategy_.ChooseHop(s, t, random);
    if (e < 0) {
      ++status.null_events;  // rejected by the strategy (thinning)
      schedule(carrier);
      continue;
    }
    if (e < graph_.EdgeOffset(s) || e >= graph_.EdgeOffset(s + 1)) {
      throw std::runtime_error(
          "KMCEngine: the waiting-time strategy chose an edge of another "
          "site");
    }
    const ToFGraph::Edge& edge = graph_.getEdge(e);
    bool sink = graph_.getSite(edge.to).role == ToFGraph::Role::Sink;
    if (!sink && occupant[edge.to] >= 0) {
      ++status.null_events;  // blocked: time has advanced, nothing happens
      schedule(carrier);
      continue;
    }
    occupant[s] = -1;
    ++status.hops;
    observer.Hop(carrier, s, edge.to, edge.dR, t);
    if (sink) {
      site_of[carrier] = -1;
      --alive;
      ++status.absorbed;
      observer.Absorb(carrier, edge.to, t);
    } else {
      occupant[edge.to] = carrier;
      site_of[carrier] = edge.to;
      schedule(carrier);
    }
  }
  status.end_time = t;
  status.events = events;
  return status;
}

}  // namespace xtp
}  // namespace votca
