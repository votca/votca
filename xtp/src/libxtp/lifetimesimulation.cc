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
#include <stdexcept>

// Third party includes
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

// Local VOTCA includes
#include "votca/xtp/lifetimesimulation.h"
#include "votca/xtp/scaledsolve.h"

namespace votca {
namespace xtp {

namespace {

using DecaySums = LifetimeSimulation::DecaySums;

class LifetimeObserver : public KMCEngine::Observer {
 public:
  LifetimeObserver(const ToFGraph& graph, Index run, bool keep,
                   Eigen::VectorXd& occupation, Eigen::VectorXd& decays,
                   std::vector<LifetimeSimulation::Record>& records,
                   Index energy_hops,
                   std::vector<LifetimeSimulation::EnergySample>* energy)
      : graph_(graph),
        run_(run),
        keep_(keep),
        occupation_(occupation),
        decays_(decays),
        records_(records),
        energy_hops_(energy_hops),
        energy_(energy) {}

  void Inject(Index carrier, Index site, double t) override {
    if (carrier >= Index(site_.size())) {
      Index n = carrier + 1;
      site_.resize(n, -1);
      entry_.resize(n, 0.0);
      injected_.resize(n, 0.0);
      hops_.resize(n, 0);
      displacement_.resize(n, Eigen::Vector3d::Zero());
    }
    site_[carrier] = site;
    entry_[carrier] = t;
    injected_[carrier] = t;
  }

  void Hop(Index carrier, Index from, Index to, const Eigen::Vector3d& dR,
           double t) override {
    occupation_[from] += t - entry_[carrier];
    entry_[carrier] = t;
    site_[carrier] = to;
    displacement_[carrier] += dR;
    ++hops_[carrier];
    ++total_hops_;
    if (energy_ != nullptr && energy_hops_ > 0 &&
        total_hops_ % energy_hops_ == 0) {
      double e = 0.0;
      Index n = 0;
      for (Index s : site_) {
        if (s >= 0) {
          e += graph_.getSite(s).energy;
          ++n;
        }
      }
      energy_->push_back({t, total_hops_, n > 0 ? e / double(n) : 0.0});
    }
  }

  void Decay(Index carrier, Index site, double t) override {
    occupation_[site] += t - entry_[carrier];
    decays_[site] += 1.0;
    double lifetime = t - injected_[carrier];
    const Eigen::Vector3d& R = displacement_[carrier];
    double r2 = R.squaredNorm();
    sums_.decays++;
    sums_.lifetime += lifetime;
    sums_.lifetime2 += lifetime * lifetime;
    sums_.distance += std::sqrt(r2);
    sums_.distance2 += r2;
    sums_.distance4 += r2 * r2;
    sums_.components += R.cwiseAbs2();
    sums_.hops += double(hops_[carrier]);
    if (keep_) {
      LifetimeSimulation::Record r;
      r.run = run_;
      r.carrier = carrier;
      r.time = t;
      r.lifetime = lifetime;
      r.hops = hops_[carrier];
      r.site = site;
      r.displacement = R;
      records_.push_back(r);
    }
    site_[carrier] = -1;
  }

  void Finish(double end) {
    for (Index c = 0; c < Index(site_.size()); ++c) {
      if (site_[c] >= 0) {
        occupation_[site_[c]] += end - entry_[c];
      }
    }
  }

  const DecaySums& Sums() const { return sums_; }

 private:
  const ToFGraph& graph_;
  Index run_;
  bool keep_;
  Eigen::VectorXd& occupation_;
  Eigen::VectorXd& decays_;
  std::vector<LifetimeSimulation::Record>& records_;
  Index energy_hops_;
  std::vector<LifetimeSimulation::EnergySample>* energy_;
  std::vector<Index> site_;
  std::vector<double> entry_;
  std::vector<double> injected_;
  std::vector<Index> hops_;
  std::vector<Eigen::Vector3d> displacement_;
  Index total_hops_ = 0;
  DecaySums sums_;
};

}  // namespace

LifetimeSimulation::LifetimeSimulation(const ToFGraph& graph,
                                       std::vector<double> decay_rates,
                                       const Settings& settings)
    : graph_(graph),
      decay_rates_(std::move(decay_rates)),
      settings_(settings),
      own_rates_(std::make_unique<StaticRates>(graph)),
      own_strategy_(std::make_unique<ExponentialWaitingTime>()) {
  Setup(*own_strategy_, *own_rates_);
}

LifetimeSimulation::LifetimeSimulation(const ToFGraph& graph,
                                       std::vector<double> decay_rates,
                                       const Settings& settings,
                                       const TimeDependentRates& rates,
                                       WaitingTimeStrategy& strategy)
    : graph_(graph), decay_rates_(std::move(decay_rates)), settings_(settings) {
  if (&rates.Graph() != &graph_) {
    throw std::runtime_error(
        "LifetimeSimulation: the rates belong to a different graph");
  }
  Setup(strategy, rates);
}

void LifetimeSimulation::Setup(WaitingTimeStrategy& strategy,
                               const TimeDependentRates& rates) {
  Check();
  strategy.Prepare(rates);
  engine_ = std::make_unique<KMCEngine>(graph_, strategy);
  KMCEngine::DecaySettings decay;
  decay.rates = &decay_rates_;
  decay.replace = true;
  engine_->setDecay(decay);
}

void LifetimeSimulation::Check() const {
  if (Index(decay_rates_.size()) != graph_.NumberOfSites()) {
    throw std::runtime_error(
        "LifetimeSimulation: need one decay rate per site");
  }
  for (Index i = 0; i < graph_.NumberOfSites(); ++i) {
    double k = decay_rates_[i];
    if (!graph_.getSite(i).ignored && !(k > 0.0 && std::isfinite(k))) {
      throw std::runtime_error("LifetimeSimulation: the decay rate of site " +
                               std::to_string(i) +
                               " must be positive and finite");
    }
  }
  if (settings_.runs < 1 || settings_.decays < settings_.runs) {
    throw std::runtime_error(
        "LifetimeSimulation: need at least one run and one decay per run");
  }
  if (settings_.carriers < 1 ||
      settings_.carriers > Index(graph_.SourceSites().size())) {
    throw std::runtime_error(
        "LifetimeSimulation: the number of carriers must be between 1 and "
        "the number of injection sites (" +
        std::to_string(graph_.SourceSites().size()) + ")");
  }
  if (settings_.energy_hops < 0) {
    throw std::runtime_error("LifetimeSimulation: energy_hops < 0");
  }
}

LifetimeSimulation::RunData LifetimeSimulation::SingleRun(
    Index run, Index decays, std::vector<EnergySample>* energy) const {
  const Index n = graph_.NumberOfSites();
  RunData data;
  data.occupation = Eigen::VectorXd::Zero(n);
  data.decays = Eigen::VectorXd::Zero(n);
  LifetimeObserver observer(graph_, run, settings_.keep_records,
                            data.occupation, data.decays, data.records,
                            settings_.energy_hops, energy);
  KMCEngine::Random random(settings_.seed, std::uint64_t(run));
  KMCEngine::Limits limits;
  limits.max_decays = decays;
  limits.max_wall_time = settings_.max_wall_time;
  std::vector<double> injection(settings_.carriers, 0.0);
  data.status = engine_->Run(injection, random, observer, limits);
  observer.Finish(data.status.end_time);
  data.sums = observer.Sums();
  return data;
}

LifetimeSimulation::Result LifetimeSimulation::Simulate() const {
  const Index runs = settings_.runs;
  const Index n = graph_.NumberOfSites();
  Result result;
  result.decay_fraction = Eigen::VectorXd::Zero(n);
  result.occupation = Eigen::VectorXd::Zero(n);
  DecaySums total;

  // blocks of runs in parallel; reduction in run order (thread independent)
  const Index block = 64;
  std::vector<RunData> data(block);
  for (Index start = 0; start < runs; start += block) {
    Index count = std::min(block, runs - start);
#pragma omp parallel for schedule(dynamic)
    for (Index k = 0; k < count; ++k) {
      Index run = start + k;
      Index decays =
          settings_.decays / runs + (run < settings_.decays % runs ? 1 : 0);
      data[k] = SingleRun(run, decays, run == 0 ? &result.energy : nullptr);
    }
    for (Index k = 0; k < count; ++k) {
      const RunData& d = data[k];
      result.occupation += d.occupation;
      result.decay_fraction += d.decays;
      result.simulated_time += d.status.end_time;
      result.hop_events += d.status.hops;
      result.null_events += d.status.null_events;
      if (d.status.reached_wall_time) {
        ++result.wall_time_stops;
      }
      const DecaySums& s = d.sums;
      total.decays += s.decays;
      total.lifetime += s.lifetime;
      total.lifetime2 += s.lifetime2;
      total.distance += s.distance;
      total.distance2 += s.distance2;
      total.distance4 += s.distance4;
      total.components += s.components;
      total.hops += s.hops;
      result.records.insert(result.records.end(), d.records.begin(),
                            d.records.end());
    }
  }
  if (total.decays == 0) {
    throw std::runtime_error("LifetimeSimulation: no carrier decayed");
  }
  double N = double(total.decays);
  result.decays = total.decays;
  auto error = [N](double sum, double sum2) {
    if (N < 2.0) {
      return 0.0;
    }
    double mean = sum / N;
    double var = std::max(sum2 / N - mean * mean, 0.0) * N / (N - 1.0);
    return std::sqrt(var / N);
  };
  result.lifetime = total.lifetime / N;
  result.lifetime_error = error(total.lifetime, total.lifetime2);
  result.mean_distance = total.distance / N;
  result.mean_distance_error = error(total.distance, total.distance2);
  result.squared_displacement = total.distance2 / N;
  result.squared_displacement_error = error(total.distance2, total.distance4);
  result.squared_components = total.components / N;
  result.hops = total.hops / N;
  result.decay_fraction /= N;
  if (result.simulated_time > 0.0) {
    result.occupation /= result.simulated_time;
  }
  return result;
}

LifetimeSimulation::ExactResult LifetimeSimulation::Exact(
    const ToFGraph& graph, const std::vector<double>& decay_rates,
    Index direct_limit, double tolerance) {
  const Index n = graph.NumberOfSites();
  if (Index(decay_rates.size()) != n) {
    throw std::runtime_error("LifetimeSimulation::Exact: one rate per site");
  }
  using SpMat = Eigen::SparseMatrix<double>;
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(graph.NumberOfEdges() + n);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(n, 3);
  for (Index i = 0; i < n; ++i) {
    double diag = decay_rates[i] + graph.EscapeRate(i);
    if (!(diag > 0.0)) {
      if (graph.getSite(i).ignored) {
        diag = 1.0;  // never visited
      } else {
        throw std::runtime_error("LifetimeSimulation::Exact: site " +
                                 std::to_string(i) +
                                 " neither decays nor hops");
      }
    }
    triplets.emplace_back(i, i, diag);
    for (const ToFGraph::Edge& e : graph.OutEdges(i)) {
      triplets.emplace_back(i, e.to, -e.rate);
      b.row(i) += e.rate * e.dR.transpose();
    }
  }
  SpMat A(n, n);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();
  SpMat At = A.transpose();
  At.makeCompressed();

  ExactResult result;
  result.converged = true;
  auto solve = [&](const SpMat& M, const Eigen::MatrixXd& rhs) {
    Eigen::MatrixXd x(rhs.rows(), rhs.cols());
    if (n <= direct_limit) {
      result.solver = "sparse LU";
      Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
      lu.compute(M);
      if (lu.info() != Eigen::Success) {
        throw std::runtime_error(
            "LifetimeSimulation::Exact: LU factorization failed");
      }
      x = lu.solve(rhs);
    } else {
      result.solver = "BiCGSTAB with ILUT";
      Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double>> bicg;
      bicg.setTolerance(tolerance);
      bicg.setMaxIterations(100000);
      bicg.compute(M);  // M outlives bicg
      for (Index c = 0; c < rhs.cols(); ++c) {
        // scaled right-hand side (see ScaledSolve); a zero one gives zero
        bool solved = false;
        x.col(c) = ScaledSolve(bicg, Eigen::VectorXd(rhs.col(c)), &solved);
        if (solved && bicg.iterations() >= bicg.maxIterations()) {
          result.converged = false;
        }
      }
    }
    for (Index c = 0; c < rhs.cols(); ++c) {
      double scale = (M * x.col(c)).cwiseAbs().maxCoeff() +
                     rhs.col(c).cwiseAbs().maxCoeff();
      double err = (M * x.col(c) - rhs.col(c)).cwiseAbs().maxCoeff() /
                   std::max(scale, std::numeric_limits<double>::min());
      result.backward_error = std::max(result.backward_error, err);
    }
    return x;
  };

  // first moments: lifetime and mean displacement
  Eigen::MatrixXd rhs1(n, 4);
  rhs1.col(0).setOnes();
  rhs1.rightCols(3) = b;
  Eigen::MatrixXd x1 = solve(A, rhs1);
  Eigen::VectorXd tau = x1.col(0);
  Eigen::MatrixXd m = x1.rightCols(3);
  // second moments per component
  Eigen::MatrixXd rhs2 = Eigen::MatrixXd::Zero(n, 3);
  for (Index i = 0; i < n; ++i) {
    for (const ToFGraph::Edge& e : graph.OutEdges(i)) {
      rhs2.row(i) +=
          e.rate * (e.dR.cwiseAbs2().transpose() +
                    2.0 * e.dR.transpose().cwiseProduct(m.row(e.to)));
    }
  }
  Eigen::MatrixXd s2 = solve(A, rhs2);
  // start distribution: uniform on the injection sites
  const std::vector<Index>& sources = graph.SourceSites();
  Eigen::VectorXd p0 = Eigen::VectorXd::Zero(n);
  for (Index s : sources) {
    p0[s] = 1.0 / double(sources.size());
  }
  Eigen::VectorXd g = solve(At, p0);

  result.lifetime = p0.dot(tau);
  result.squared_components = (p0.transpose() * s2).transpose();
  result.squared_displacement = result.squared_components.sum();
  Eigen::VectorXd k = Eigen::Map<const Eigen::VectorXd>(decay_rates.data(), n);
  result.decay_probability = k.cwiseProduct(g);
  result.occupation = g / g.sum();
  result.converged = result.converged && result.backward_error < 1e-10;
  return result;
}

}  // namespace xtp
}  // namespace votca
