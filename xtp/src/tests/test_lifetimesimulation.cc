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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE lifetimesimulation_test

// Standard includes
#include <cmath>
#include <cstdint>
#include <unordered_map>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/eigen.h"
#include "votca/xtp/kmcengine.h"
#include "votca/xtp/latticemodel.h"
#include "votca/xtp/lifetimesimulation.h"

using namespace votca::xtp;
using votca::Index;

namespace {
const double nm = votca::tools::conv::nm2bohr;
const double ev = votca::tools::conv::ev2hrt;
const double kT = votca::tools::conv::kB * 300.0 * ev;

// 3x3x3 periodic cubic lattice, nearest neighbours
void BuildLattice(Topology& top, LatticeModel::Disorder disorder) {
  LatticeModel::Options opt;
  opt.dims = {3, 3, 3};
  opt.spacing = 1.0 * nm;
  opt.cutoff = 1.1 * nm;
  opt.J2 = 1.79e-4 * ev * ev;
  opt.decay_length = 1.0 * nm;
  opt.reorganization = 0.23 * ev;
  opt.mean_energy = -0.76 * ev;
  opt.sigma = 0.1 * ev;
  opt.disorder = disorder;
  opt.seed = 1;
  opt.states = {QMStateType::Singlet};
  LatticeModel(opt).BuildTopology(top);
}

std::vector<double> DecayRates(Index n) {
  std::vector<double> k(n);
  for (Index i = 0; i < n; ++i) {
    k[i] = 1.0 / (1e-10 * (1.0 + 0.6 * std::sin(1.7 * double(i))));
  }
  return k;
}

ToFGraph Bulk(const Topology& top) {
  return ToFGraph(top, ToFGraph::BulkOptions("*", ""), QMStateType::Singlet, kT,
                  Eigen::Vector3d::Zero());
}

// N carriers with exclusion, every decay replaced on a uniformly chosen free
// site: mean lifetime N / (stationary decay rate) by Little's law, from all
// configurations (n <= 32)
double ExclusionLifetime(const ToFGraph& g, const std::vector<double>& k,
                         Index N) {
  const Index n = g.NumberOfSites();
  std::vector<std::uint32_t> configs;
  for (std::uint64_t c = 0; c < (std::uint64_t(1) << n); ++c) {
    if (Index(__builtin_popcountll(c)) == N) {
      configs.push_back(std::uint32_t(c));
    }
  }
  std::unordered_map<std::uint32_t, Index> index;
  for (Index q = 0; q < Index(configs.size()); ++q) {
    index[configs[q]] = q;
  }
  const Index M = Index(configs.size());
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(M, M);
  Eigen::VectorXd decay = Eigen::VectorXd::Zero(M);
  for (Index q = 0; q < M; ++q) {
    std::uint32_t c = configs[q];
    for (Index s = 0; s < n; ++s) {
      if (!(c & (1u << s))) {
        continue;
      }
      for (const ToFGraph::Edge& e : g.OutEdges(s)) {
        if (!(c & (1u << e.to))) {
          L(q, index.at((c & ~(1u << s)) | (1u << e.to))) += e.rate;
        }
      }
      decay[q] += k[s];
      std::uint32_t without = c & ~(1u << s);
      double free = double(n - N + 1);
      for (Index j = 0; j < n; ++j) {
        if (!(without & (1u << j))) {
          L(q, index.at(without | (1u << j))) += k[s] / free;
        }
      }
    }
    L(q, q) -= L.row(q).sum();
  }
  Eigen::MatrixXd A = L.transpose();
  A.row(0).setOnes();
  Eigen::VectorXd b = Eigen::VectorXd::Zero(M);
  b(0) = 1.0;
  Eigen::VectorXd p = A.partialPivLu().solve(b);
  return double(N) / p.dot(decay);
}

LifetimeSimulation::Settings Settings(Index carriers, Index runs,
                                      Index decays) {
  LifetimeSimulation::Settings s;
  s.seed = 5;
  s.runs = runs;
  s.carriers = carriers;
  s.decays = decays;
  return s;
}

void CheckValue(double kmc, double error, double exact, double max_rel) {
  BOOST_CHECK_GT(error, 0.0);
  BOOST_CHECK_LT(std::abs(kmc - exact), 4.5 * error);
  BOOST_CHECK_LT(error, max_rel * std::abs(exact));
}
}  // namespace

BOOST_AUTO_TEST_SUITE(lifetimesimulation_test)

BOOST_AUTO_TEST_CASE(exact_uniform_decay_on_ordered_lattice) {
  // equal energies and decay rate k: lifetime 1/k, and a random walk with
  // total hop rate W makes on average W/k hops of length a, uncorrelated
  Topology top;
  BuildLattice(top, LatticeModel::Disorder::None);
  ToFGraph graph = Bulk(top);
  const double k = 1e10;
  std::vector<double> rates(graph.NumberOfSites(), k);
  LifetimeSimulation::ExactResult r = LifetimeSimulation::Exact(graph, rates);
  double W = graph.EscapeRate(0);
  BOOST_CHECK(r.converged);
  BOOST_CHECK_CLOSE(r.lifetime, 1.0 / k, 1e-9);
  BOOST_CHECK_CLOSE(r.squared_displacement, W / k * nm * nm, 1e-8);
  for (Index d = 0; d < 3; ++d) {
    BOOST_CHECK_CLOSE(r.squared_components[d], W / (3.0 * k) * nm * nm, 1e-8);
  }
  BOOST_CHECK_CLOSE(r.decay_probability.sum(), 1.0, 1e-10);
  BOOST_CHECK_SMALL(r.decay_probability.maxCoeff() - 1.0 / 27.0, 1e-12);
  // the iterative solver gives the same
  LifetimeSimulation::ExactResult it =
      LifetimeSimulation::Exact(graph, rates, 0);
  BOOST_CHECK(it.converged);
  BOOST_CHECK_CLOSE(it.lifetime, r.lifetime, 1e-8);
  BOOST_CHECK_CLOSE(it.squared_displacement, r.squared_displacement, 1e-8);
}

BOOST_AUTO_TEST_CASE(one_carrier_against_exact) {
  Topology top;
  BuildLattice(top, LatticeModel::Disorder::Uncorrelated);
  ToFGraph graph = Bulk(top);
  std::vector<double> k = DecayRates(graph.NumberOfSites());
  LifetimeSimulation::ExactResult exact = LifetimeSimulation::Exact(graph, k);
  BOOST_CHECK(exact.converged);
  LifetimeSimulation sim(graph, k, Settings(1, 4, 200000));
  LifetimeSimulation::Result r = sim.Simulate();
  BOOST_CHECK_EQUAL(r.decays, 200000);
  CheckValue(r.lifetime, r.lifetime_error, exact.lifetime, 0.01);
  CheckValue(r.squared_displacement, r.squared_displacement_error,
             exact.squared_displacement, 0.01);
  BOOST_CHECK_LT(
      (r.decay_fraction - exact.decay_probability).cwiseAbs().maxCoeff(),
      0.005);
  BOOST_CHECK_LT((r.occupation - exact.occupation).cwiseAbs().maxCoeff(),
                 0.005);
  BOOST_CHECK_CLOSE(r.occupation.sum(), 1.0, 1e-8);
  // the diffusion length is at least the mean distance
  BOOST_CHECK_GE(std::sqrt(r.squared_displacement), r.mean_distance);
}

BOOST_AUTO_TEST_CASE(exclusion_against_exact) {
  // the previous kmclifetime gave 10 % too short lifetimes with 3 excitons
  // on this lattice (time advanced only by the full escape rates)
  Topology top;
  BuildLattice(top, LatticeModel::Disorder::Uncorrelated);
  ToFGraph graph = Bulk(top);
  std::vector<double> k = DecayRates(graph.NumberOfSites());
  for (Index N : {2, 25}) {
    double exact = ExclusionLifetime(graph, k, N);
    // nearly full lattice: most attempts are blocked, fewer decays suffice
    Index decays = (N == 2) ? 200000 : 30000;
    LifetimeSimulation sim(graph, k, Settings(N, 4, decays));
    LifetimeSimulation::Result r = sim.Simulate();
    BOOST_TEST_MESSAGE("N " << N << " exact " << exact << " kmc " << r.lifetime
                            << " +- " << r.lifetime_error << " hops "
                            << r.hop_events << " null " << r.null_events);
    CheckValue(r.lifetime, r.lifetime_error, exact, 0.01);
    BOOST_CHECK_CLOSE(r.occupation.sum(), double(N), 1e-6);
    BOOST_CHECK_GT(r.null_events, 0);
  }
}

BOOST_AUTO_TEST_CASE(records_energy_and_reproducibility) {
  Topology top;
  BuildLattice(top, LatticeModel::Disorder::Uncorrelated);
  ToFGraph graph = Bulk(top);
  std::vector<double> k = DecayRates(graph.NumberOfSites());
  LifetimeSimulation::Settings s = Settings(3, 3, 1000);
  s.keep_records = true;
  s.energy_hops = 10;
  LifetimeSimulation sim(graph, k, s);
  LifetimeSimulation::Result a = sim.Simulate();
  LifetimeSimulation::Result b = sim.Simulate();
  BOOST_CHECK_EQUAL(a.hop_events, b.hop_events);
  BOOST_CHECK_EQUAL(a.lifetime, b.lifetime);
  BOOST_REQUIRE_EQUAL(a.records.size(), 1000);
  // decays split 334/333/333 over the runs, in run order
  BOOST_CHECK_EQUAL(a.records[333].run, 0);
  BOOST_CHECK_EQUAL(a.records[334].run, 1);
  double sum = 0.0;
  for (const LifetimeSimulation::Record& r : a.records) {
    BOOST_CHECK_GT(r.lifetime, 0.0);
    BOOST_CHECK_LE(r.lifetime, r.time);
    sum += r.lifetime;
  }
  BOOST_CHECK_CLOSE(sum / 1000.0, a.lifetime, 1e-10);
  BOOST_CHECK(!a.energy.empty());
  BOOST_CHECK_EQUAL(a.energy.front().hops, 10);
}

BOOST_AUTO_TEST_CASE(engine_decay_without_replacement) {
  Topology top;
  BuildLattice(top, LatticeModel::Disorder::Uncorrelated);
  ToFGraph graph = Bulk(top);
  std::vector<double> k = DecayRates(graph.NumberOfSites());
  StaticRates rates(graph);
  ExponentialWaitingTime strategy;
  strategy.Prepare(rates);
  KMCEngine engine(graph, strategy);
  KMCEngine::DecaySettings decay;
  decay.rates = &k;
  engine.setDecay(decay);
  KMCEngine::Random random(1, 0);
  KMCEngine::Observer observer;
  KMCEngine::RunStatus st = engine.Run(std::vector<double>(5, 0.0), random,
                                       observer, KMCEngine::Limits());
  BOOST_CHECK_EQUAL(st.injected, 5);
  BOOST_CHECK_EQUAL(st.decayed, 5);
  BOOST_CHECK(!st.censored);
  std::vector<double> wrong(3, 1.0);
  decay.rates = &wrong;
  BOOST_CHECK_THROW(engine.setDecay(decay), std::runtime_error);
  decay.rates = nullptr;
  decay.replace = true;
  BOOST_CHECK_THROW(engine.setDecay(decay), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(invalid_setups) {
  Topology top;
  BuildLattice(top, LatticeModel::Disorder::None);
  ToFGraph graph = Bulk(top);
  std::vector<double> k(graph.NumberOfSites(), 1e9);
  k[4] = 0.0;
  BOOST_CHECK_THROW(LifetimeSimulation(graph, k, Settings(1, 1, 10)),
                    std::runtime_error);
  k[4] = 1e9;
  BOOST_CHECK_THROW(LifetimeSimulation(graph, k, Settings(28, 1, 10)),
                    std::runtime_error);
  BOOST_CHECK_THROW(LifetimeSimulation(graph, k, Settings(1, 5, 4)),
                    std::runtime_error);
  k.pop_back();
  BOOST_CHECK_THROW(LifetimeSimulation(graph, k, Settings(1, 1, 10)),
                    std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
