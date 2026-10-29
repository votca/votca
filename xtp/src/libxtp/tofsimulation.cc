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
#include <cmath>
#include <functional>
#include <numeric>
#include <stdexcept>

// Local VOTCA includes
#include "votca/xtp/eigen.h"
#include "votca/xtp/tofsimulation.h"

namespace votca {
namespace xtp {

namespace {
class RunObserver : public KMCEngine::Observer {
 public:
  RunObserver(const Eigen::Vector3d& normal, Index carriers,
              std::function<void(double, double)> displacement,
              std::function<void(double)> arrival)
      : normal_(normal),
        inject_time_(carriers, 0.0),
        displacement_(std::move(displacement)),
        arrival_(std::move(arrival)) {}
  void Inject(Index carrier, Index, double t) override {
    inject_time_[carrier] = t;
  }
  void Hop(Index, Index, Index, const Eigen::Vector3d& dR, double t) override {
    displacement_(t, dR.dot(normal_));
  }
  void Absorb(Index carrier, Index, double t) override {
    if (absorbed_ == 0) {
      first_ = t;
    }
    ++absorbed_;
    transit_sum_ += t - inject_time_[carrier];
    arrival_(t);
  }
  double First() const { return first_; }
  double MeanTransit() const {
    return absorbed_ > 0 ? transit_sum_ / double(absorbed_) : 0.0;
  }

 private:
  Eigen::Vector3d normal_;
  std::vector<double> inject_time_;
  std::function<void(double, double)> displacement_;
  std::function<void(double)> arrival_;
  Index absorbed_ = 0;
  double first_ = 0.0;
  double transit_sum_ = 0.0;
};
}  // namespace

ToFSimulation::ToFSimulation(const ToFGraph& graph, const Settings& settings)
    : graph_(graph),
      settings_(settings),
      own_rates_(std::make_unique<StaticRates>(graph)),
      own_strategy_(std::make_unique<ExponentialWaitingTime>()) {
  own_strategy_->Prepare(*own_rates_);
  engine_ = std::make_unique<KMCEngine>(graph_, *own_strategy_);
  strategy_name_ = own_strategy_->Name();
  Setup();
}

ToFSimulation::ToFSimulation(const ToFGraph& graph, const Settings& settings,
                             const TimeDependentRates& rates,
                             WaitingTimeStrategy& strategy)
    : graph_(graph), settings_(settings) {
  if (&rates.Graph() != &graph_) {
    throw std::runtime_error(
        "ToFSimulation: the rates belong to a different graph");
  }
  strategy.Prepare(rates);
  engine_ = std::make_unique<KMCEngine>(graph_, strategy);
  strategy_name_ = strategy.Name();
  Setup();
}

void ToFSimulation::Setup() {
  if (graph_.isBulk()) {
    throw std::runtime_error(
        "ToFSimulation: needs a ToF graph with sinks, not a bulk graph without "
        "sinks");
  }

  if (settings_.runs < 1) {
    throw std::runtime_error("ToFSimulation: runs must be positive");
  }
  if (!(settings_.tmin > 0.0) || !(settings_.tmax > settings_.tmin) ||
      settings_.bins_per_decade < 1) {
    throw std::runtime_error(
        "ToFSimulation: need 0 < tmin < tmax and bins_per_decade >= 1");
  }
  Index bins = Index(std::ceil(double(settings_.bins_per_decade) *
                               std::log10(settings_.tmax / settings_.tmin)));
  bins = std::max<Index>(bins, 1);
  edges_.resize(bins + 1);
  for (Index b = 0; b <= bins; ++b) {
    edges_[b] = settings_.tmin *
                std::pow(10.0, double(b) / double(settings_.bins_per_decade));
  }
}

Index ToFSimulation::Bin(double t) const {
  if (t < edges_.front()) {
    return -1;
  }
  if (t >= edges_.back()) {
    return Index(edges_.size()) - 1;  // overflow
  }
  return Index(std::upper_bound(edges_.begin(), edges_.end(), t) -
               edges_.begin()) -
         1;
}

ToFSimulation::RunData ToFSimulation::SingleRun(Index run,
                                                std::vector<double>& current,
                                                std::vector<double>& arrivals,
                                                double& before,
                                                double& after) const {
  // two streams per run: the injection times and the dynamics
  KMCEngine::Random pulse_random(settings_.seed, 2 * std::uint64_t(run));
  KMCEngine::Random random(settings_.seed, 2 * std::uint64_t(run) + 1);
  std::vector<double> times =
      KMCEngine::InjectionTimes(settings_.pulses, pulse_random);
  const Index nbins = Index(edges_.size()) - 1;
  RunObserver observer(
      graph_.Normal(), Index(times.size()),
      [&](double t, double dx) {
        Index b = Bin(t);
        if (b < 0) {
          before += dx;
        } else if (b >= nbins) {
          after += dx;
        } else {
          current[b] += dx;
        }
      },
      [&](double t) {
        Index b = Bin(t);
        if (b < nbins) {
          arrivals[std::max<Index>(b, 0)] += 1.0;
        }
      });
  RunData data;
  data.status = engine_->Run(times, random, observer, settings_.limits);
  data.first_arrival = observer.First();
  data.mean_transit = observer.MeanTransit();
  return data;
}

ToFSimulation::Statistic ToFSimulation::Analyze(std::vector<double> samples) {
  Statistic s;
  s.samples = Index(samples.size());
  if (samples.empty()) {
    return s;
  }
  double n = double(samples.size());
  s.mean = std::accumulate(samples.begin(), samples.end(), 0.0) / n;
  double var = std::accumulate(
      samples.begin(), samples.end(), 0.0,
      [&s](double acc, double x) { return acc + (x - s.mean) * (x - s.mean); });
  var = (samples.size() > 1) ? var / (n - 1.0) : 0.0;
  s.standard_deviation = std::sqrt(var);
  s.standard_error = s.standard_deviation / std::sqrt(n);
  std::sort(samples.begin(), samples.end(), std::greater<double>());
  Index top = std::max<Index>(1, Index(std::ceil(0.01 * n)));
  double top_sum = std::accumulate(samples.begin(), samples.begin() + top, 0.0);
  s.top_percent_share = (s.mean > 0.0) ? top_sum / (s.mean * n) : 0.0;
  return s;
}

ToFSimulation::Result ToFSimulation::Simulate() const {
  const Index nbins = Index(edges_.size()) - 1;
  const Index runs = settings_.runs;
  Result result;
  result.runs.resize(runs);
  result.bin_edges = edges_;
  std::vector<double> current(nbins, 0.0);
  std::vector<double> arrivals(nbins, 0.0);
  double before = 0.0;
  double after = 0.0;

  // blocks of runs in parallel; reduction in run order (thread independent)
  const Index block = 256;
  std::vector<std::vector<double>> run_current(block);
  std::vector<std::vector<double>> run_arrivals(block);
  std::vector<double> run_before(block);
  std::vector<double> run_after(block);
  for (Index start = 0; start < runs; start += block) {
    Index count = std::min(block, runs - start);
#pragma omp parallel for schedule(dynamic)
    for (Index k = 0; k < count; ++k) {
      run_current[k].assign(nbins, 0.0);
      run_arrivals[k].assign(nbins, 0.0);
      run_before[k] = 0.0;
      run_after[k] = 0.0;
      result.runs[start + k] =
          SingleRun(start + k, run_current[k], run_arrivals[k], run_before[k],
                    run_after[k]);
    }
    for (Index k = 0; k < count; ++k) {
      for (Index b = 0; b < nbins; ++b) {
        current[b] += run_current[k][b];
        arrivals[b] += run_arrivals[k][b];
      }
      before += run_before[k];
      after += run_after[k];
    }
  }

  std::vector<double> first;
  std::vector<double> transit;
  for (const RunData& r : result.runs) {
    result.injected += r.status.injected;
    result.blocked_injections += r.status.blocked_injections;
    result.hops += r.status.hops;
    result.null_events += r.status.null_events;
    if (r.status.censored) {
      ++result.censored;
      continue;
    }
    if (r.status.absorbed > 0) {
      first.push_back(r.first_arrival);
      transit.push_back(r.mean_transit);
    }
  }
  result.first = Analyze(first);
  result.transit = Analyze(transit);

  double norm = graph_.Length() * double(std::max<Index>(result.injected, 1));
  result.current.resize(nbins);
  result.collected.resize(nbins);
  double cumulative = 0.0;
  for (Index b = 0; b < nbins; ++b) {
    result.current[b] = current[b] / (norm * (edges_[b + 1] - edges_[b]));
    cumulative += arrivals[b];
    result.collected[b] =
        cumulative / double(std::max<Index>(result.injected, 1));
  }
  result.displacement_before = before / norm;
  result.displacement_after = after / norm;
  return result;
}

}  // namespace xtp
}  // namespace votca
