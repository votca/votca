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
#include <limits>
#include <stdexcept>

// Third party includes
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

// Local VOTCA includes
#include "votca/xtp/bulksimulation.h"
#include "votca/xtp/scaledsolve.h"

namespace votca {
namespace xtp {

namespace {

/// Follows the carriers of one run: unwrapped positions, occupation times,
/// positions at the block boundaries t0 + j tau (tau doubled whenever more
/// than 2 * blocks + 1 boundaries are stored) and optional samples.
class BulkObserver : public KMCEngine::Observer {
 public:
  BulkObserver(const ToFGraph& graph, Index carriers, double start,
               Index blocks, double tau0, Eigen::VectorXd& occupation,
               double sample_interval, Index sample_hops,
               std::vector<BulkSimulation::Sample>* samples)
      : graph_(graph),
        position_(carriers, Eigen::Vector3d::Zero()),
        site_(carriers, -1),
        entry_(carriers, 0.0),
        start_(start),
        blocks_(blocks),
        tau_(tau0),
        occupation_(occupation),
        sample_interval_(sample_interval),
        sample_hops_(sample_hops),
        samples_(samples) {}

  void Inject(Index carrier, Index site, double t) override {
    Advance(t);
    position_[carrier] = graph_.getSite(site).pos;
    site_[carrier] = site;
    entry_[carrier] = t;
  }

  void Hop(Index carrier, Index from, Index to, const Eigen::Vector3d& dR,
           double t) override {
    Advance(t);
    Occupy(from, entry_[carrier], t);
    position_[carrier] += dR;
    site_[carrier] = to;
    entry_[carrier] = t;
    ++hops_;
    if (samples_ != nullptr && sample_hops_ > 0 && hops_ % sample_hops_ == 0) {
      Record(t);
    }
  }

  /// close the run at time end
  void Finish(double end) {
    Advance(end);
    for (Index c = 0; c < Index(site_.size()); ++c) {
      if (site_[c] >= 0) {
        Occupy(site_[c], entry_[c], end);
      }
    }
  }

  const std::vector<std::vector<Eigen::Vector3d>>& Boundaries() const {
    return boundaries_;
  }
  double Tau() const { return tau_; }
  const std::vector<Eigen::Vector3d>& Positions() const { return position_; }

 private:
  const ToFGraph& graph_;
  std::vector<Eigen::Vector3d> position_;
  std::vector<Index> site_;
  std::vector<double> entry_;
  double start_;
  Index blocks_;
  double tau_;
  Eigen::VectorXd& occupation_;
  double sample_interval_;
  Index sample_hops_;
  std::vector<BulkSimulation::Sample>* samples_;
  std::vector<std::vector<Eigen::Vector3d>> boundaries_;
  Index next_sample_ = 0;
  Index hops_ = 0;

  void Occupy(Index site, double from, double to) {
    from = std::max(from, start_);
    if (to > from) {
      occupation_[site] += to - from;
    }
  }

  /// record the state before an event at time t (the state is constant
  /// since the previous event)
  void Advance(double t) {
    while (start_ + double(boundaries_.size()) * tau_ < t) {
      if (Index(boundaries_.size()) == 2 * blocks_ + 1) {
        // keep the boundaries at even multiples of tau
        std::vector<std::vector<Eigen::Vector3d>> kept;
        kept.reserve(blocks_ + 1);
        for (Index j = 0; j < Index(boundaries_.size()); j += 2) {
          kept.push_back(std::move(boundaries_[j]));
        }
        boundaries_ = std::move(kept);
        tau_ *= 2.0;
        continue;
      }
      boundaries_.push_back(position_);
    }
    if (samples_ == nullptr || !(sample_interval_ > 0.0)) {
      return;
    }
    while (double(next_sample_) * sample_interval_ < t) {
      Record(double(next_sample_) * sample_interval_);
      ++next_sample_;
    }
  }

  void Record(double t) {
    BulkSimulation::Sample s;
    s.time = t;
    s.hops = hops_;
    s.positions = position_;
    double energy = 0.0;
    Index n = 0;
    for (Index site : site_) {
      if (site >= 0) {
        energy += graph_.getSite(site).energy;
        ++n;
      }
    }
    s.energy = (n > 0) ? energy / double(n) : 0.0;
    samples_->push_back(std::move(s));
  }
};

}  // namespace

BulkSimulation::BulkSimulation(const ToFGraph& graph, const Settings& settings)
    : graph_(graph),
      settings_(settings),
      own_rates_(std::make_unique<StaticRates>(graph)),
      own_strategy_(std::make_unique<ExponentialWaitingTime>()) {
  Check();
  own_strategy_->Prepare(*own_rates_);
  engine_ = std::make_unique<KMCEngine>(graph_, *own_strategy_);
}

BulkSimulation::BulkSimulation(const ToFGraph& graph, const Settings& settings,
                               const TimeDependentRates& rates,
                               WaitingTimeStrategy& strategy)
    : graph_(graph), settings_(settings) {
  Check();
  if (&rates.Graph() != &graph_) {
    throw std::runtime_error(
        "BulkSimulation: the rates belong to a different graph");
  }
  strategy.Prepare(rates);
  engine_ = std::make_unique<KMCEngine>(graph_, strategy);
}

void BulkSimulation::Check() const {
  if (!graph_.isBulk()) {
    throw std::runtime_error(
        "BulkSimulation: needs a bulk graph (ToFGraph::BulkOptions), not a ToF "
        "graph with sinks");
  }
  if (settings_.runs < 1 || settings_.blocks < 2) {
    throw std::runtime_error(
        "BulkSimulation: need at least one run and two blocks");
  }
  if (settings_.carriers < 1 ||
      settings_.carriers > Index(graph_.SourceSites().size())) {
    throw std::runtime_error(
        "BulkSimulation: the number of carriers must be between 1 and the "
        "number of injection sites (" +
        std::to_string(graph_.SourceSites().size()) + ")");
  }
  if (settings_.equilibration < 0.0 || settings_.sample_interval < 0.0 ||
      settings_.sample_hops < 0) {
    throw std::runtime_error(
        "BulkSimulation: equilibration and sampling interval must not be "
        "negative");
  }
  const KMCEngine::Limits& l = settings_.limits;
  if (l.max_time == std::numeric_limits<double>::infinity() &&
      l.max_events == std::numeric_limits<Index>::max()) {
    throw std::runtime_error(
        "BulkSimulation: a limit on the simulated time or on the number of "
        "events is required");
  }
  if (settings_.equilibration >= l.max_time) {
    throw std::runtime_error(
        "BulkSimulation: the equilibration time is not shorter than the run");
  }
}

BulkSimulation::RunData BulkSimulation::SingleRun(
    Index run, Eigen::VectorXd& occupation,
    std::vector<Sample>* trajectory) const {
  const Index n = settings_.carriers;
  const KMCEngine::Limits& limits = settings_.limits;
  double max_escape = 0.0;
  for (Index i = 0; i < graph_.NumberOfSites(); ++i) {
    max_escape = std::max(max_escape, graph_.EscapeRate(i));
  }
  // any small start works, tau is doubled as needed
  double tau0 = (max_escape > 0.0) ? 1e-3 / max_escape : 1e-18;
  BulkObserver observer(graph_, n, settings_.equilibration, settings_.blocks,
                        tau0, occupation, settings_.sample_interval,
                        settings_.sample_hops, trajectory);
  KMCEngine::Random random(settings_.seed, std::uint64_t(run));
  std::vector<double> injection(n, 0.0);

  RunData data;
  data.status = engine_->Run(injection, random, observer, limits);
  // the state is constant after the last event: if the run stopped at the
  // time limit, or because nobody could move, the window extends to it
  const KMCEngine::RunStatus& st = data.status;
  bool stopped_early = st.reached_wall_time || st.events >= limits.max_events;
  double end = (!stopped_early && std::isfinite(limits.max_time))
                   ? limits.max_time
                   : st.end_time;
  observer.Finish(end);
  data.start = settings_.equilibration;
  data.end = end;
  double window = end - data.start;
  if (!(window > 0.0)) {
    occupation.setZero();
    return data;  // the run ended before the equilibration time
  }
  occupation /= window;

  const auto& b = observer.Boundaries();
  data.displacement.resize(n);
  for (Index c = 0; c < n; ++c) {
    data.displacement[c] = observer.Positions()[c] - b.front()[c];
    data.velocity += data.displacement[c];
  }
  data.velocity /= double(n) * window;

  Index blocks = Index(b.size()) - 1;
  data.blocks = blocks;
  data.block_length = observer.Tau();
  if (blocks < 2) {
    return data;
  }
  const double tau = observer.Tau();
  Eigen::Vector3d mean = Eigen::Vector3d::Zero();
  for (Index k = 0; k < blocks; ++k) {
    for (Index c = 0; c < n; ++c) {
      mean += b[k + 1][c] - b[k][c];
    }
  }
  double samples = double(blocks * n);
  mean /= samples;
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d block_cov = Eigen::Matrix3d::Zero();
  for (Index k = 0; k < blocks; ++k) {
    Eigen::Vector3d per_block = Eigen::Vector3d::Zero();
    for (Index c = 0; c < n; ++c) {
      Eigen::Vector3d d = b[k + 1][c] - b[k][c] - mean;
      cov += d * d.transpose();
      per_block += d;
    }
    per_block /= double(n) * tau;  // velocity of the block minus the mean
    block_cov += per_block * per_block.transpose();
  }
  if (samples > 1.0) {
    data.diffusion = cov / ((samples - 1.0) * 2.0 * tau);
  }
  data.velocity_covariance = block_cov / (double(blocks - 1) * double(blocks));
  return data;
}

BulkSimulation::Result BulkSimulation::Simulate() const {
  const Index runs = settings_.runs;
  const Index sites = graph_.NumberOfSites();
  Result result;
  result.runs.resize(runs);
  result.occupation = Eigen::VectorXd::Zero(sites);

  // blocks of runs in parallel; reduction in run order (thread independent)
  const Index block = 64;
  std::vector<Eigen::VectorXd> run_occupation(block);
  for (Index start = 0; start < runs; start += block) {
    Index count = std::min(block, runs - start);
#pragma omp parallel for schedule(dynamic)
    for (Index k = 0; k < count; ++k) {
      run_occupation[k] = Eigen::VectorXd::Zero(sites);
      Index run = start + k;
      result.runs[run] = SingleRun(run, run_occupation[k],
                                   run == 0 ? &result.trajectory : nullptr);
    }
    for (Index k = 0; k < count; ++k) {
      result.occupation += run_occupation[k];
    }
  }
  result.occupation /= double(runs);

  Index valid = 0;
  for (const RunData& r : result.runs) {
    result.hops += r.status.hops;
    result.null_events += r.status.null_events;
    if (r.status.reached_wall_time) {
      ++result.wall_time_stops;
    }
    if (r.end > r.start) {
      result.velocity += r.velocity;
      result.diffusion += r.diffusion;
      ++valid;
    }
  }
  if (valid == 0) {
    throw std::runtime_error(
        "BulkSimulation: no run lasted longer than the equilibration time");
  }
  result.velocity /= double(valid);
  result.diffusion /= double(valid);
  if (valid == 1) {
    for (const RunData& r : result.runs) {
      if (r.end > r.start) {
        result.velocity_error =
            r.velocity_covariance.diagonal().cwiseMax(0.0).cwiseSqrt();
      }
    }
    return result;
  }
  Eigen::Vector3d var_v = Eigen::Vector3d::Zero();
  Eigen::Matrix3d var_d = Eigen::Matrix3d::Zero();
  for (const RunData& r : result.runs) {
    if (r.end > r.start) {
      var_v += (r.velocity - result.velocity).cwiseAbs2();
      var_d += (r.diffusion - result.diffusion).cwiseAbs2();
    }
  }
  double norm = double(valid) * double(valid - 1);
  result.velocity_error = (var_v / norm).cwiseSqrt();
  result.diffusion_error = (var_d / norm).cwiseSqrt();
  return result;
}

namespace {
/// strongly connected components (iterative Tarjan) over edges with
/// positive rate; returns the component of every site
std::vector<Index> Components(const ToFGraph& g, Index& count) {
  const Index n = g.NumberOfSites();
  std::vector<Index> index(n, -1), low(n, 0), comp(n, -1);
  std::vector<char> on_stack(n, 0);
  std::vector<Index> stack;
  std::vector<std::pair<Index, Index>> call;  // (site, next edge)
  Index counter = 0;
  count = 0;
  for (Index root = 0; root < n; ++root) {
    if (index[root] >= 0) {
      continue;
    }
    call.emplace_back(root, g.EdgeOffset(root));
    index[root] = low[root] = counter++;
    stack.push_back(root);
    on_stack[root] = 1;
    while (!call.empty()) {
      Index v = call.back().first;
      Index& e = call.back().second;
      if (e < g.EdgeOffset(v + 1)) {
        const ToFGraph::Edge& edge = g.getEdge(e++);
        if (!(edge.rate > 0.0)) {
          continue;
        }
        Index w = edge.to;
        if (index[w] < 0) {
          index[w] = low[w] = counter++;
          stack.push_back(w);
          on_stack[w] = 1;
          call.emplace_back(w, g.EdgeOffset(w));
        } else if (on_stack[w]) {
          low[v] = std::min(low[v], index[w]);
        }
        continue;
      }
      if (low[v] == index[v]) {
        Index w;
        do {
          w = stack.back();
          stack.pop_back();
          on_stack[w] = 0;
          comp[w] = count;
        } while (w != v);
        ++count;
      }
      call.pop_back();
      if (!call.empty()) {
        Index u = call.back().first;
        low[u] = std::min(low[u], low[v]);
      }
    }
  }
  return comp;
}
}  // namespace

BulkSimulation::ExactResult BulkSimulation::Exact(const ToFGraph& graph,
                                                  Index direct_limit,
                                                  double tolerance) {
  ExactResult result;
  const Index n = graph.NumberOfSites();
  result.occupation = Eigen::VectorXd::Zero(n);
  Index ncomp = 0;
  std::vector<Index> comp = Components(graph, ncomp);
  // closed classes: no edge with positive rate leaves the component; ignored
  // sites are isolated and never visited, so they do not count
  std::vector<char> closed(ncomp, 1);
  std::vector<char> relevant(ncomp, 0);
  for (Index i = 0; i < n; ++i) {
    if (!graph.getSite(i).ignored) {
      relevant[comp[i]] = 1;
    }
    for (const ToFGraph::Edge& e : graph.OutEdges(i)) {
      if (e.rate > 0.0 && comp[e.to] != comp[i]) {
        closed[comp[i]] = 0;
      }
    }
  }
  std::vector<Index> classes;
  for (Index c = 0; c < ncomp; ++c) {
    if (closed[c] && relevant[c]) {
      classes.push_back(c);
    }
  }
  if (classes.size() != 1) {
    result.message = std::to_string(classes.size()) +
                     " closed classes of sites (e.g. isolated sites or "
                     "clusters without connection): the stationary "
                     "distribution depends on the start";
    return result;
  }
  // local numbering of the closed class
  std::vector<Index> local(n, -1);
  std::vector<Index> sites;
  for (Index i = 0; i < n; ++i) {
    if (comp[i] == classes[0]) {
      local[i] = Index(sites.size());
      sites.push_back(i);
    }
  }
  const Index m = Index(sites.size());
  result.sites = m;
  if (m == 1) {
    result.available = true;
    result.converged = true;
    result.solver = "trivial";
    result.occupation[sites[0]] = 1.0;
    return result;
  }

  using SpMat = Eigen::SparseMatrix<double>;
  // generator rows (from): L_ii = -W_i, L_ij = w_ij
  std::vector<Eigen::Triplet<double>> tl;
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(m, 3);
  for (Index a = 0; a < m; ++a) {
    Index i = sites[a];
    double out = 0.0;
    for (const ToFGraph::Edge& e : graph.OutEdges(i)) {
      if (e.rate > 0.0) {
        tl.emplace_back(a, local[e.to], e.rate);
        out += e.rate;
        b.row(a) += e.rate * e.dR.transpose();
      }
    }
    tl.emplace_back(a, a, -out);
  }
  SpMat L(m, m);
  L.setFromTriplets(tl.begin(), tl.end());
  double scale = 0.0;
  for (Index a = 0; a < m; ++a) {
    scale = std::max(scale, -L.coeff(a, a));
  }

  result.converged = true;
  auto solve = [&](SpMat M, const Eigen::MatrixXd& rhs) {
    M.makeCompressed();
    Eigen::MatrixXd x(rhs.rows(), rhs.cols());
    if (m <= direct_limit) {
      result.solver = "sparse LU";
      Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
      lu.compute(M);
      if (lu.info() != Eigen::Success) {
        throw std::runtime_error("BulkSimulation::Exact: LU failed");
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
      Eigen::VectorXd r = M * x.col(c) - rhs.col(c);
      double denom = (M * x.col(c)).cwiseAbs().maxCoeff() +
                     rhs.col(c).cwiseAbs().maxCoeff();
      if (denom > 0.0) {
        result.backward_error =
            std::max(result.backward_error, r.cwiseAbs().maxCoeff() / denom);
      }
    }
    return x;
  };

  // One equation of L^T pi = 0 and of L chi = v - b is dependent; it is
  // replaced by fixing the value at one site ("pin"), which keeps the
  // matrices sparse (a normalization or gauge row would be dense and fill
  // the LU and ILUT factors). The pin is the site of lowest energy, where pi
  // is large, so that the other values are not huge. Afterwards pi is
  // normalized and chi is shifted to the gauge pi^T chi = 0 (L 1 = 0).
  Index pin = 0;
  for (Index a = 1; a < m; ++a) {
    if (graph.getSite(sites[a]).energy < graph.getSite(sites[pin]).energy) {
      pin = a;
    }
  }
  auto pinned = [pin](SpMat M) {
    M.prune([pin](Index row, Index, double) { return row != pin; });
    M.coeffRef(pin, pin) = 1.0;
    return M;
  };
  SpMat A = pinned(SpMat(L.transpose()) / scale);
  Eigen::MatrixXd e0 = Eigen::MatrixXd::Zero(m, 1);
  e0(pin, 0) = 1.0;
  Eigen::VectorXd pi = solve(A, e0).col(0);
  pi /= pi.sum();
  Eigen::Vector3d v = (pi.transpose() * b).transpose();
  SpMat B = pinned(L / scale);
  Eigen::MatrixXd rhs = (Eigen::VectorXd::Ones(m) * v.transpose() - b) / scale;
  rhs.row(pin).setZero();
  Eigen::MatrixXd chi = solve(B, rhs);
  for (Index c = 0; c < 3; ++c) {
    chi.col(c).array() -= pi.dot(chi.col(c));
  }

  Eigen::Matrix3d D = Eigen::Matrix3d::Zero();
  for (Index a = 0; a < m; ++a) {
    for (const ToFGraph::Edge& e : graph.OutEdges(sites[a])) {
      if (e.rate > 0.0) {
        Eigen::Vector3d g =
            e.dR + (chi.row(local[e.to]) - chi.row(a)).transpose();
        D += 0.5 * pi(a) * e.rate * g * g.transpose();
      }
    }
  }
  for (Index a = 0; a < m; ++a) {
    result.occupation[sites[a]] = pi(a);
  }
  result.velocity = v;
  result.diffusion = D;
  result.available = true;
  result.converged = result.converged && result.backward_error < 1e-10;
  return result;
}

}  // namespace xtp
}  // namespace votca
