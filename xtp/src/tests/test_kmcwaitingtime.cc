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

#define BOOST_TEST_MODULE kmcwaitingtime_test

// Standard includes
#include <cmath>
#include <vector>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/graphrandomwalk.h"
#include "votca/xtp/kmcrates.h"
#include "votca/xtp/kmcwaitingtime.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/tofsimulation.h"

using namespace votca::xtp;
using votca::Index;

namespace {
const double nm = votca::tools::conv::nm2bohr;
const double ev = votca::tools::conv::ev2hrt;
const double kT = votca::tools::conv::kB * 300.0 * ev;
const double Pi = votca::tools::conv::Pi;

/*
 * Time-dependent test rates: all static rates of the graph multiplied by
 * one periodic factor s(t) = factor (1 + a cos(2 pi t / T + phi)).
 * Uniform scaling only reparametrizes time: with
 *   S(t) = int_0^t s = factor (t + a T/(2 pi) (sin(2 pi t/T + phi) - sin phi))
 * a trajectory of the static chain reaching an event at static time theta
 * reaches it at S^-1(theta) under s(t). This gives exact references for any
 * time-dependent strategy, also with several carriers and exclusion.
 */
class ModulatedRates final : public TimeDependentRates {
 public:
  ModulatedRates(const ToFGraph& graph, double factor, double a, double T,
                 double phi)
      : TimeDependentRates(graph), factor_(factor), a_(a), T_(T), phi_(phi) {}
  double Scale(double t) const {
    return factor_ * (1.0 + a_ * std::cos(2.0 * Pi * t / T_ + phi_));
  }
  double Rate(Index edge, double t) const override {
    return graph_.getEdge(edge).rate * Scale(t);
  }
  double UpperBound(Index edge) const override {
    return graph_.getEdge(edge).rate * factor_ * (1.0 + std::abs(a_));
  }
  double Period() const override { return T_; }
  // static time theta -> time under the modulation
  double MapStaticTime(double theta) const {
    auto S = [this](double t) {
      return factor_ *
             (t + a_ * T_ / (2.0 * Pi) *
                      (std::sin(2.0 * Pi * t / T_ + phi_) - std::sin(phi_)));
    };
    double lo = 0.0;
    double hi = theta / (factor_ * (1.0 - std::abs(a_))) + T_;
    for (Index it = 0; it < 200; ++it) {
      double mid = 0.5 * (lo + hi);
      (S(mid) < theta ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
  }

 private:
  double factor_;
  double a_;
  double T_;
  double phi_;
};

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

struct Sample {
  double mean = 0.0;
  double standard_error = 0.0;
};

// mean and standard error of the first arrivals of uncensored runs, mapped
Sample FirstArrivals(const ToFSimulation::Result& r,
                     const ModulatedRates* map = nullptr) {
  std::vector<double> x;
  for (const ToFSimulation::RunData& d : r.runs) {
    if (!d.status.censored && d.status.absorbed > 0) {
      x.push_back(map ? map->MapStaticTime(d.first_arrival) : d.first_arrival);
    }
  }
  ToFSimulation::Statistic s = ToFSimulation::Analyze(x);
  return {s.mean, s.standard_error};
}

void CheckAgree(const Sample& a, const Sample& b) {
  double se = std::sqrt(a.standard_error * a.standard_error +
                        b.standard_error * b.standard_error);
  BOOST_CHECK_LT(std::abs(a.mean - b.mean), 4.5 * se);
  BOOST_CHECK_LT(se, 0.02 * std::abs(b.mean));
}

ToFSimulation::Settings TwoCarrierSettings(Index runs, std::uint64_t seed) {
  ToFSimulation::Settings settings;
  settings.runs = runs;
  settings.seed = seed;
  settings.pulses[0].carriers = 2;
  return settings;
}

// exact reference for modulated rates: static (exponential) runs mapped
// through the inverse integrated modulation
Sample MappedStaticReference(const ToFGraph& graph, const ModulatedRates& rates,
                             Index runs) {
  StaticRates fixed(graph);
  ExponentialWaitingTime exponential;
  ToFSimulation reference(graph, TwoCarrierSettings(runs, 99), fixed,
                          exponential);
  return FirstArrivals(reference.Simulate(), &rates);
}
}  // namespace

BOOST_AUTO_TEST_SUITE(kmcwaitingtime_test)

BOOST_AUTO_TEST_CASE(exponential_needs_static_rates) {
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  ModulatedRates rates(graph, 1.0, 0.5, 1e-10, 0.0);
  ExponentialWaitingTime exponential;
  BOOST_CHECK_THROW(exponential.Prepare(rates), std::runtime_error);
  BOOST_CHECK(MakeWaitingTimeStrategy("thinning")->Name() == "thinning");
  BOOST_CHECK_THROW(MakeWaitingTimeStrategy("bisection"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(thinning_static_rates_against_grw) {
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  GraphRandomWalk grw(graph);
  double first = grw.TwoCarriers(GraphRandomWalk::Observable::First).tau_mean;
  double transit =
      grw.TwoCarriers(GraphRandomWalk::Observable::Transit).tau_mean;
  StaticRates rates(graph);
  ThinningWaitingTime thinning;
  ToFSimulation sim(graph, TwoCarrierSettings(30000, 3), rates, thinning);
  ToFSimulation::Result r = sim.Simulate();
  BOOST_CHECK_LT(std::abs(r.first.mean - first), 4.5 * r.first.standard_error);
  BOOST_CHECK_LT(std::abs(r.transit.mean - transit),
                 4.5 * r.transit.standard_error);
}

BOOST_AUTO_TEST_CASE(thinning_constant_factor) {
  // rates halved, bounds at the full rates: half of the proposals are
  // rejected and all times double
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  double first = GraphRandomWalk(graph)
                     .TwoCarriers(GraphRandomWalk::Observable::First)
                     .tau_mean;
  ModulatedRates half(graph, 0.5, 0.0, 1.0, 0.0);
  ThinningWaitingTime thinning;
  ToFSimulation sim(graph, TwoCarrierSettings(30000, 4), half, thinning);
  ToFSimulation::Result r = sim.Simulate();
  BOOST_CHECK_LT(std::abs(r.first.mean - 2.0 * first),
                 4.5 * r.first.standard_error);
}

BOOST_AUTO_TEST_CASE(thinning_modulated_rates_against_time_map) {
  // two carriers with exclusion, rates modulated by 1 + 0.8 cos(2 pi t / T
  // + phi) with T a few first-arrival times: starting at the maximum (phi =
  // 0) or the minimum (phi = pi) of the rates changes the times strongly
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  double first = GraphRandomWalk(graph)
                     .TwoCarriers(GraphRandomWalk::Observable::First)
                     .tau_mean;
  for (double phi : {0.0, Pi}) {
    ModulatedRates rates(graph, 1.0, 0.8, 4.0 * first, phi);
    Sample reference = MappedStaticReference(graph, rates, 30000);
    // the modulation matters: the mapped mean differs from the static one
    BOOST_CHECK_GT(std::abs(reference.mean - first),
                   10.0 * reference.standard_error);
    ThinningWaitingTime thinning;
    ToFSimulation sim(graph, TwoCarrierSettings(30000, 5), rates, thinning);
    CheckAgree(FirstArrivals(sim.Simulate()), reference);
  }
}

BOOST_AUTO_TEST_CASE(inversion_is_a_template) {
  // remove this test case when IntegratedRateInversion is implemented and
  // enable the ones below
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  ModulatedRates rates(graph, 1.0, 0.5, 1e-10, 0.0);
  IntegratedRateInversion inversion{IntegratedRateInversion::Settings()};
  BOOST_CHECK_THROW(inversion.Prepare(rates), std::logic_error);
  // static rates have no period
  StaticRates fixed(graph);
  BOOST_CHECK_THROW(inversion.Prepare(fixed), std::runtime_error);
}

// Acceptance tests for IntegratedRateInversion (BEP nhKMC): enable when
// implemented (remove the disabled() decorator or run with
// --run_test=kmcwaitingtime_test/inversion_*).

BOOST_AUTO_TEST_CASE(inversion_constant_rates_against_grw,
                     *boost::unit_test::disabled()) {
  // modulation amplitude 0: the inversion must reproduce the static chain
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  GraphRandomWalk grw(graph);
  double first = grw.TwoCarriers(GraphRandomWalk::Observable::First).tau_mean;
  double transit =
      grw.TwoCarriers(GraphRandomWalk::Observable::Transit).tau_mean;
  ModulatedRates rates(graph, 1.0, 0.0, 0.3 * first, 0.0);
  IntegratedRateInversion inversion{IntegratedRateInversion::Settings()};
  ToFSimulation sim(graph, TwoCarrierSettings(30000, 6), rates, inversion);
  ToFSimulation::Result r = sim.Simulate();
  BOOST_CHECK_LT(std::abs(r.first.mean - first), 4.5 * r.first.standard_error);
  BOOST_CHECK_LT(std::abs(r.transit.mean - transit),
                 4.5 * r.transit.standard_error);
}

BOOST_AUTO_TEST_CASE(inversion_modulated_rates_against_time_map,
                     *boost::unit_test::disabled()) {
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  double first = GraphRandomWalk(graph)
                     .TwoCarriers(GraphRandomWalk::Observable::First)
                     .tau_mean;
  // periods from much shorter to longer than the first-arrival time
  for (double T : {0.01 * first, first, 10.0 * first}) {
    ModulatedRates rates(graph, 1.0, 0.8, T, 1.0);
    Sample reference = MappedStaticReference(graph, rates, 30000);
    IntegratedRateInversion inversion{IntegratedRateInversion::Settings()};
    ToFSimulation sim(graph, TwoCarrierSettings(30000, 7), rates, inversion);
    CheckAgree(FirstArrivals(sim.Simulate()), reference);
  }
}

BOOST_AUTO_TEST_CASE(inversion_against_thinning,
                     *boost::unit_test::disabled()) {
  // general periodic rates (here: modulated); thinning is the reference for
  // rates that are not a uniform scaling, e.g. Marcus rates in an
  // oscillating field
  Topology top;
  BuildFiveSites(top);
  ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  ModulatedRates rates(graph, 1.0, 0.6, 2e-10, 0.5);
  ThinningWaitingTime thinning;
  IntegratedRateInversion inversion{IntegratedRateInversion::Settings()};
  ToFSimulation a(graph, TwoCarrierSettings(30000, 8), rates, thinning);
  ToFSimulation b(graph, TwoCarrierSettings(30000, 9), rates, inversion);
  CheckAgree(FirstArrivals(b.Simulate()), FirstArrivals(a.Simulate()));
}

BOOST_AUTO_TEST_SUITE_END()
