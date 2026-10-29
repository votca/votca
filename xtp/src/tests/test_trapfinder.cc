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

#define BOOST_TEST_MODULE trapfinder_test

// Standard includes
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/stochasticnetwork.h"
#include "votca/xtp/stochasticnetworkfit.h"
#include "votca/xtp/tofgraph.h"
#include "votca/xtp/trapfinder.h"

using namespace votca::xtp;
using votca::Index;

namespace {
const double nm = votca::tools::conv::nm2bohr;
const double ev = votca::tools::conv::ev2hrt;
const double kT300 = 300.0 * votca::tools::conv::kB * ev;

/// a small amorphous network for holes with weak disorder
StochasticNetwork::Options SmallNetwork(Index sites, std::uint64_t seed) {
  StochasticNetwork::Options opt;
  opt.sites = sites;
  double l = std::cbrt(double(sites) / 1.3);
  opt.box = Eigen::Vector3d(l, l, l) * nm;
  opt.points.hard_core = 0.37 * nm;
  opt.points.radius_mean = 0.15 * nm;
  opt.points.radius_variance = 0.008 * nm * nm;
  opt.connection_r = {0.5 * nm, 1.0 * nm, 1.5 * nm};
  opt.connection_p = {1.0, 1.0, 0.0};
  StochasticNetwork::Carrier h;
  h.state = QMStateType(QMStateType::Hole);
  h.energy.sigma = 0.05 * ev;
  h.coupling.mean_offset = -0.3;
  h.coupling.mean_slope = -5.0;
  h.coupling.var_offset = 0.3;
  h.U_nX_nN = 0.15 * ev;
  h.U_xN_xX = 0.15 * ev;
  opt.carriers = {h};
  opt.seed = seed;
  return opt;
}

TrapFinder::Result FindTraps(const StochasticNetwork& net) {
  Topology top;
  net.BuildTopology(top);
  ToFGraph graph(top, ToFGraph::BulkOptions("", ""),
                 QMStateType(QMStateType::Hole), kT300,
                 Eigen::Vector3d::Zero());
  return TrapFinder::Find(graph);
}
}  // namespace

BOOST_AUTO_TEST_SUITE(trapfinder_test)

BOOST_AUTO_TEST_CASE(no_traps_in_weak_gaussian_disorder) {
  for (std::uint64_t seed = 1; seed <= 6; ++seed) {
    StochasticNetwork net(SmallNetwork(300, seed));
    TrapFinder::Result r = FindTraps(net);
    BOOST_CHECK(r.traps.empty());
    BOOST_REQUIRE(r.eigenvalues.size() >= 3);
    BOOST_CHECK(r.eigenvalues[1] > 0.0);
    double sum = 0.0;
    for (double p : r.occupation) {
      BOOST_CHECK(p >= 0.0);
      sum += p;
    }
    BOOST_CHECK_CLOSE(sum, 1.0, 1e-8);
  }
}

BOOST_AUTO_TEST_CASE(finds_inserted_single_site_and_cluster) {
  StochasticNetwork::Options opt = SmallNetwork(300, 3);
  StochasticNetwork::Trap single;
  single.density = 1.0 / 300.0;
  single.offsets = {{-0.6 * ev}};
  StochasticNetwork::Trap cluster;
  cluster.density = 1.0 / 300.0;
  cluster.offsets = {{-0.45 * ev}, {-0.45 * ev}, {-0.45 * ev}};
  opt.traps = {single};
  {
    StochasticNetwork net(opt);
    BOOST_REQUIRE_EQUAL(net.Traps().size(), 1);
    TrapFinder::Result r = FindTraps(net);
    BOOST_REQUIRE_EQUAL(r.traps.size(), 1);
    BOOST_CHECK(r.traps[0].release_rate * 10.0 < r.eigenvalues[2]);
    BOOST_CHECK(r.traps[0].sites == net.Traps()[0].sites);
    BOOST_CHECK(r.traps[0].weight > 0.99);
    BOOST_CHECK(r.traps[0].release_rate < r.traps[0].capture_rate);
    BOOST_CHECK(r.traps[0].release_rate < r.traps[0].hop_rate);
  }
  opt.traps = {cluster};
  {
    StochasticNetwork net(opt);
    BOOST_REQUIRE_EQUAL(net.Traps().size(), 1);
    std::vector<Index> inserted = net.Traps()[0].sites;
    BOOST_REQUIRE_EQUAL(inserted.size(), 3);
    TrapFinder::Result r = FindTraps(net);
    BOOST_REQUIRE_EQUAL(r.traps.size(), 1);
    std::vector<Index> found = r.traps[0].sites;
    std::sort(found.begin(), found.end());
    std::sort(inserted.begin(), inserted.end());
    BOOST_CHECK(found == inserted);
  }
}

BOOST_AUTO_TEST_CASE(insertion_counts_and_factor) {
  StochasticNetwork::Options opt = SmallNetwork(200, 5);
  StochasticNetwork::Trap t;
  t.density = 2.5 / 200.0;  // 2 or 3 traps
  t.offsets = {{-0.5 * ev}, {-0.3 * ev}};
  opt.traps = {t};
  StochasticNetwork net(opt);
  Index count = Index(net.Traps().size());
  BOOST_CHECK(count == 2 || count == 3);
  std::vector<char> used(200, 0);
  for (const auto& tr : net.Traps()) {
    BOOST_REQUIRE_EQUAL(tr.sites.size(), 2);
    BOOST_CHECK_CLOSE(net.SiteEnergies(0)[tr.sites[0]], -0.5 * ev, 1e-10);
    BOOST_CHECK_CLOSE(net.SiteEnergies(0)[tr.sites[1]], -0.3 * ev, 1e-10);
    for (Index i : tr.sites) {
      BOOST_CHECK(!used[i]);  // traps do not overlap
      used[i] = 1;
    }
  }
  // the traps change neither the positions nor the pairs
  StochasticNetwork::Options none = opt;
  none.trap_factor = 0.0;
  StochasticNetwork plain(none);
  BOOST_CHECK(plain.Traps().empty());
  BOOST_REQUIRE_EQUAL(plain.Pairs().size(), net.Pairs().size());
  BOOST_CHECK_EQUAL(plain.Pairs()[7].J2[0], net.Pairs()[7].J2[0]);
  BOOST_CHECK(plain.Positions()[11] == net.Positions()[11]);
  // ten times the density
  StochasticNetwork::Options more = opt;
  more.trap_factor = 10.0;
  StochasticNetwork many(more);
  BOOST_CHECK(Index(many.Traps().size()) >= 25);
}

BOOST_AUTO_TEST_CASE(fit_turns_traps_into_templates) {
  StochasticNetwork::Options opt = SmallNetwork(800, 7);
  StochasticNetwork::Trap t;
  t.density = 1.0 / 800.0;
  t.offsets = {{-0.6 * ev}};
  opt.traps = {t};
  StochasticNetwork net(opt);
  std::vector<QMStateType> carriers = {QMStateType(QMStateType::Hole)};
  StochasticNetworkFit::Data data =
      StochasticNetworkFit::FromNetwork(net, carriers, opt.box);
  StochasticNetworkFit::Settings s;
  s.carriers = carriers;
  s.fit_points = false;
  s.points = opt.points;
  s.l_b = {2, 4};
  s.l_c = {20, 40};
  TrapFinder::Result r = FindTraps(net);
  BOOST_REQUIRE_EQUAL(r.traps.size(), 1);
  s.traps.push_back({0, r.traps[0].sites, r.traps[0].release_rate});
  StochasticNetworkFit::Result res = StochasticNetworkFit(data, s).Fit();
  BOOST_REQUIRE_EQUAL(res.options.traps.size(), 1);
  const auto& tr = res.options.traps[0];
  BOOST_CHECK_CLOSE(tr.density, 1.0 / 800.0, 1e-10);
  BOOST_REQUIRE_EQUAL(tr.offsets.size(), 1);
  // relative to the bulk mean, which is close to the model mean 0
  BOOST_CHECK_SMALL(tr.offsets[0][0] / ev - (-0.6), 0.02);
  // the bulk sigma is not inflated by the trap
  BOOST_CHECK_CLOSE(res.energies[0].sigma, 0.05 * ev, 8.0);
  // without excluding the trap the sample sigma is larger
  StochasticNetworkFit::Settings s0 = s;
  s0.traps.clear();
  StochasticNetworkFit::Result res0 = StochasticNetworkFit(data, s0).Fit();
  BOOST_CHECK(res0.energies[0].sigma > res.energies[0].sigma * 1.05);
}

BOOST_AUTO_TEST_CASE(traps_in_parameter_file) {
  StochasticNetwork::Options opt = SmallNetwork(100, 1);
  StochasticNetwork::Carrier e = opt.carriers[0];
  e.state = QMStateType(QMStateType::Electron);
  opt.carriers.push_back(e);
  StochasticNetwork::Trap t;
  t.carrier = 1;
  t.density = 0.002;
  t.release_rate = 79.0;
  t.offsets = {{0.1 * ev, -0.7 * ev}, {0.0, -0.4 * ev}};
  opt.traps = {t};
  const std::string file = "trapfinder_test.xml";
  {
    std::ofstream out(file);
    out << StochasticNetwork::ToXML(opt);
  }
  votca::tools::Property prop;
  prop.LoadFromXML(file);
  std::remove(file.c_str());
  StochasticNetwork::Options back =
      StochasticNetwork::FromProperty(prop.get("stochasticnetwork"));
  BOOST_REQUIRE_EQUAL(back.traps.size(), 1);
  const auto& b = back.traps[0];
  BOOST_CHECK_EQUAL(b.carrier, 1);
  BOOST_CHECK_CLOSE(b.density, 0.002, 1e-8);
  BOOST_CHECK_CLOSE(b.release_rate, 79.0, 1e-8);
  BOOST_REQUIRE_EQUAL(b.offsets.size(), 2);
  BOOST_REQUIRE_EQUAL(b.offsets[1].size(), 2);
  BOOST_CHECK_CLOSE(b.offsets[0][1], -0.7 * ev, 1e-6);
  BOOST_CHECK_SMALL(b.offsets[1][0], 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
