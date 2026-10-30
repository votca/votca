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

#define BOOST_TEST_MODULE kmcengine_test

// Standard includes
#include <cmath>
#include <numeric>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/eigen.h"
#include "votca/xtp/graphrandomwalk.h"
#include "votca/xtp/kmcengine.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/tofsimulation.h"

using namespace votca::xtp;
using votca::Index;

namespace {
const double nm = votca::tools::conv::nm2bohr;
const double ev = votca::tools::conv::ev2hrt;
const double kT = votca::tools::conv::kB * 300.0 * ev;
const double Vm = ev / (1e9 * nm);

// chain of n sites along x, source = site 0, sink = site n-1
void BuildChain(Topology& top, Index n) {
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
  box.diagonal() << double(n + 2) * nm, 5 * nm, 5 * nm;
  top.setBox(box);
  for (Index i = 0; i < n; ++i) {
    Segment& seg = top.AddSegment("s");
    seg.push_back(
        Atom(0, "C", Eigen::Vector3d(0.5 + double(i), 2.5, 2.5) * nm));
    seg.setEMpoles(QMStateType::Hole, 0.03 * std::sin(double(i)) * ev);
    seg.setU_nX_nN(0.1 * ev, QMStateType::Hole);
    seg.setU_xN_xX(0.1 * ev, QMStateType::Hole);
  }
  for (Index i = 0; i + 1 < n; ++i) {
    QMPair& pair = top.NBList().Add(top.Segments()[i], top.Segments()[i + 1],
                                    Eigen::Vector3d(1, 0, 0) * nm);
    pair.setJeff2(1e-4 * ev * ev, QMStateType::Hole);
  }
}

ToFGraph::Options ChainRegions() {
  ToFGraph::Options opt;
  opt.source_thickness = 1.0 * nm;
  opt.sink_thickness = 2.6 * nm;
  return opt;
}

// five-site graph of Fig. 2 in arXiv:2405.15836, sources 0 and 1, sink 4
void BuildFiveSites(Topology& top) {
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
  box.diagonal() << 3.5 * nm, 5 * nm, 5 * nm;
  top.setBox(box);
  std::vector<Eigen::Vector3d> pos = {{0.5, 2.5, 2.5},
                                      {0.5, 1.5, 2.5},
                                      {1.5, 2.5, 2.5},
                                      {1.5, 1.5, 2.5},
                                      {2.5, 2.0, 2.5}};
  std::vector<double> energy = {0.02, -0.05, 0.08, -0.11, 0.0};
  for (Index i = 0; i < 5; ++i) {
    Segment& seg = top.AddSegment("s");
    seg.push_back(Atom(0, "C", pos[i] * nm));
    seg.setEMpoles(QMStateType::Hole, energy[i] * ev);
    seg.setU_nX_nN(0.1 * ev, QMStateType::Hole);
    seg.setU_xN_xX(0.1 * ev, QMStateType::Hole);
  }
  std::vector<std::array<Index, 2>> edges = {{0, 1}, {0, 2}, {1, 3},
                                             {2, 3}, {2, 4}, {3, 4}};
  double J2 = 1e-4;
  for (const auto& e : edges) {
    const Segment& a = top.Segments()[e[0]];
    const Segment& b = top.Segments()[e[1]];
    QMPair& pair = top.NBList().Add(a, b, b.getPos() - a.getPos());
    pair.setJeff2(J2 * ev * ev, QMStateType::Hole);
    J2 *= 0.7;
  }
}

ToFGraph::Options FiveSiteRegions() {
  ToFGraph::Options opt;
  opt.source_thickness = 1.0 * nm;
  opt.sink_thickness = 1.2 * nm;
  return opt;
}

// |kmc - exact| within 4.5 standard errors (false alarm ~ 1e-5)
void CheckStatistic(const ToFSimulation::Statistic& s, double exact) {
  BOOST_CHECK_GT(s.samples, 0);
  BOOST_CHECK_LT(std::abs(s.mean - exact), 4.5 * s.standard_error);
  // and the standard error is small enough for this to mean something
  BOOST_CHECK_LT(s.standard_error, 0.02 * exact);
}
}  // namespace

BOOST_AUTO_TEST_SUITE(kmcengine_test)

BOOST_AUTO_TEST_CASE(random_numbers) {
  KMCEngine::Random a(7, 3);
  KMCEngine::Random b(7, 3);
  KMCEngine::Random c(7, 4);
  double sum = 0.0;
  bool all_equal = true;
  bool any_different = false;
  for (Index i = 0; i < 1000; ++i) {
    double x = a.Uniform();
    all_equal = all_equal && (x == b.Uniform());
    any_different = any_different || (x != c.Uniform());
    BOOST_CHECK(x > 0.0 && x <= 1.0);
    sum += x;
  }
  BOOST_CHECK(all_equal);
  BOOST_CHECK(any_different);
  BOOST_CHECK_CLOSE(sum / 1000.0, 0.5, 5.0);
  std::vector<Index> counts(5, 0);
  for (Index i = 0; i < 5000; ++i) {
    Index k = a.Integer(5);
    BOOST_REQUIRE(k >= 0 && k < 5);
    ++counts[k];
  }
  for (Index k : counts) {
    BOOST_CHECK(k > 850 && k < 1150);
  }
}

BOOST_AUTO_TEST_CASE(pulses) {
  KMCEngine::Random random(1, 0);
  KMCEngine::Pulse p1;
  p1.carriers = 3;
  KMCEngine::Pulse p2;
  p2.start = 2e-9;
  p2.carriers = 50;
  p2.shape = KMCEngine::Pulse::Shape::Gaussian;
  p2.width = 1e-9;
  KMCEngine::Pulse p3;
  p3.start = 1e-6;
  p3.carriers = 20;
  p3.shape = KMCEngine::Pulse::Shape::Rectangular;
  p3.width = 1e-7;
  std::vector<double> t = KMCEngine::InjectionTimes({p1, p2, p3}, random);
  BOOST_REQUIRE_EQUAL(t.size(), 73);
  BOOST_CHECK(std::is_sorted(t.begin(), t.end()));
  BOOST_CHECK_EQUAL(t[0], 0.0);
  BOOST_CHECK_EQUAL(t[2], 0.0);
  BOOST_CHECK_GE(t[3], 0.0);
  for (Index k = 53; k < 73; ++k) {
    BOOST_CHECK(t[k] >= 1e-6 && t[k] <= 1.1e-6);
  }
  BOOST_CHECK(KMCEngine::ShapeFromString("gaussian") ==
              KMCEngine::Pulse::Shape::Gaussian);
  BOOST_CHECK_THROW(KMCEngine::ShapeFromString("box"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(single_carrier_chain_against_grw) {
  Topology top;
  BuildChain(top, 10);
  for (Eigen::Vector3d field :
       {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(2e7 * Vm, 0, 0)}) {
    ToFGraph graph(top, ChainRegions(), QMStateType::Hole, kT, field);
    GraphRandomWalk grw(graph);
    GraphRandomWalk::Result exact = grw.SingleCarrier();
    ToFSimulation::Settings settings;
    settings.runs = 20000;
    settings.tmin = 1e-15;
    settings.tmax = 1e-3;
    ToFSimulation simulation(graph, settings);
    ToFSimulation::Result r = simulation.Simulate();
    BOOST_CHECK_EQUAL(r.censored, 0);
    CheckStatistic(r.first, exact.tau_mean);
    CheckStatistic(r.transit, exact.tau_mean);
    BOOST_CHECK_CLOSE(r.first.standard_deviation, exact.tau_std_injection, 5.0);
    // every carrier moves from x = 0.5 to 9.5 nm, L = 9 nm: the current
    // integrates to exactly one
    double integral = r.displacement_before + r.displacement_after;
    for (Index b = 0; b < Index(r.current.size()); ++b) {
      integral += r.current[b] * (r.bin_edges[b + 1] - r.bin_edges[b]);
    }
    BOOST_CHECK_CLOSE(integral, 1.0, 1e-8);
    BOOST_CHECK_CLOSE(r.collected.back(), 1.0, 1e-12);
  }
}

BOOST_AUTO_TEST_CASE(two_carriers_with_blocking_against_grw) {
  // five sites, two carriers: about 60 % of the proposals hit the occupied
  // site, so wrong waiting times would show immediately
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  GraphRandomWalk grw(graph);
  double first = grw.TwoCarriers(GraphRandomWalk::Observable::First).tau_mean;
  double transit =
      grw.TwoCarriers(GraphRandomWalk::Observable::Transit).tau_mean;
  ToFSimulation::Settings settings;
  settings.runs = 100000;
  settings.pulses[0].carriers = 2;
  ToFSimulation simulation(graph, settings);
  ToFSimulation::Result r = simulation.Simulate();
  BOOST_CHECK_EQUAL(r.censored, 0);
  BOOST_CHECK_EQUAL(r.injected, 200000);
  BOOST_CHECK_GT(double(r.null_events), 0.5 * double(r.hops));
  CheckStatistic(r.first, first);
  CheckStatistic(r.transit, transit);
}

BOOST_AUTO_TEST_CASE(thread_independence) {
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  ToFSimulation::Settings settings;
  settings.runs = 700;  // not a multiple of the block size
  settings.pulses[0].carriers = 2;
  ToFSimulation simulation(graph, settings);
  Index threads = OPENMP::getMaxThreads();
  OPENMP::setMaxThreads(1);
  ToFSimulation::Result one = simulation.Simulate();
  OPENMP::setMaxThreads(3);
  ToFSimulation::Result three = simulation.Simulate();
  OPENMP::setMaxThreads(threads);
  BOOST_CHECK_EQUAL(one.first.mean, three.first.mean);
  BOOST_CHECK_EQUAL(one.transit.mean, three.transit.mean);
  BOOST_CHECK_EQUAL(one.hops, three.hops);
  BOOST_CHECK(one.current == three.current);
  for (Index k = 0; k < settings.runs; ++k) {
    BOOST_CHECK_EQUAL(one.runs[k].first_arrival, three.runs[k].first_arrival);
  }
}

BOOST_AUTO_TEST_CASE(limits_and_blocked_injection) {
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  // three carriers at once on two source sites: one injection is dropped
  ToFSimulation::Settings settings;
  settings.runs = 50;
  settings.pulses[0].carriers = 3;
  ToFSimulation::Result r = ToFSimulation(graph, settings).Simulate();
  BOOST_CHECK_EQUAL(r.blocked_injections, 50);
  BOOST_CHECK_EQUAL(r.injected, 100);
  // a time limit far below the hopping time: every run is censored
  settings.pulses[0].carriers = 1;
  settings.limits.max_time = 1e-16;
  r = ToFSimulation(graph, settings).Simulate();
  BOOST_CHECK_EQUAL(r.censored, 50);
  BOOST_CHECK_EQUAL(r.transit.samples, 0);
  // an event limit
  settings.limits.max_time = std::numeric_limits<double>::infinity();
  settings.limits.max_events = 1;
  r = ToFSimulation(graph, settings).Simulate();
  BOOST_CHECK_GT(r.censored, 0);
}

BOOST_AUTO_TEST_CASE(statistics) {
  ToFSimulation::Statistic s = ToFSimulation::Analyze({1.0, 2.0, 3.0, 4.0});
  BOOST_CHECK_EQUAL(s.samples, 4);
  BOOST_CHECK_CLOSE(s.mean, 2.5, 1e-12);
  BOOST_CHECK_CLOSE(s.standard_deviation, std::sqrt(5.0 / 3.0), 1e-12);
  BOOST_CHECK_CLOSE(s.standard_error, std::sqrt(5.0 / 3.0) / 2.0, 1e-12);
  BOOST_CHECK_CLOSE(s.top_percent_share, 0.4, 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
