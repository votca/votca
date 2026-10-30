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

#define BOOST_TEST_MODULE stochasticnetwork_test

// Standard includes
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/segment.h"
#include "votca/xtp/stochasticnetwork.h"
#include "votca/xtp/stochasticnetworkfit.h"

using namespace votca::xtp;
using votca::Index;

namespace {
const double nm = votca::tools::conv::nm2bohr;
const double ev = votca::tools::conv::ev2hrt;

/// parameters in the range of the amorphous systems of the paper
StochasticNetwork::Options TestOptions(Index sites) {
  StochasticNetwork::Options opt;
  opt.sites = sites;
  double l = std::cbrt(double(sites) / 1.3);  // about 1.3 sites per nm^3
  opt.box = Eigen::Vector3d(l, l, l) * nm;
  opt.points.hard_core = 0.37 * nm;
  opt.points.radius_mean = 0.15 * nm;
  opt.points.radius_variance = 0.008 * nm * nm;
  opt.connection_r = {0.5 * nm, 1.0 * nm, 1.5 * nm, 2.0 * nm};
  opt.connection_p = {1.0, 1.0, 0.4, 0.0};
  StochasticNetwork::Carrier h;
  h.state = QMStateType(QMStateType::Hole);
  h.energy.mean = 0.1 * ev;
  h.energy.sigma = 0.09 * ev;
  h.energy.w_a = 0.3;
  h.energy.w_b = 0.3;
  h.energy.l_b = 6;
  h.energy.l_c = 60;
  h.coupling.mean_offset = -0.3;
  h.coupling.mean_slope = -5.0;
  h.coupling.var_offset = 0.45;
  h.coupling.var_slope = 1.5;
  h.U_nX_nN = 0.12 * ev;
  h.U_xN_xX = 0.18 * ev;
  StochasticNetwork::Carrier e = h;
  e.state = QMStateType(QMStateType::Electron);
  e.energy.sigma = 0.1 * ev;
  e.energy.l_b = 3;
  e.energy.l_c = 40;
  e.coupling.mean_offset = -0.6;
  opt.carriers = {h, e};
  opt.energy_correlation = -0.8;
  opt.coupling_correlation = 0.5;
  opt.seed = 3;
  return opt;
}

std::vector<Eigen::Vector3d> RandomPoints(Index n, const Eigen::Vector3d& box,
                                          std::uint64_t seed) {
  ModelRandom r(seed, 0);
  std::vector<Eigen::Vector3d> p(n);
  for (auto& x : p) {
    for (Index d = 0; d < 3; ++d) {
      x[d] = (1.0 - r.Uniform()) * box[d];
    }
  }
  return p;
}

double MinImage(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                const Eigen::Vector3d& box) {
  Eigen::Vector3d d = b - a;
  for (Index k = 0; k < 3; ++k) {
    d[k] -= box[k] * std::round(d[k] / box[k]);
  }
  return d.norm();
}
}  // namespace

BOOST_AUTO_TEST_SUITE(stochasticnetwork_test)

BOOST_AUTO_TEST_CASE(random_numbers) {
  ModelRandom a(5, 1);
  ModelRandom b(5, 1);
  ModelRandom c(5, 2);
  double sa = 0.0, sc = 0.0;
  for (Index k = 0; k < 10; ++k) {
    double x = a.Uniform();
    BOOST_CHECK_EQUAL(x, b.Uniform());
    BOOST_CHECK(x > 0.0 && x <= 1.0);
    sa += x;
    sc += c.Uniform();
  }
  BOOST_CHECK(sa != sc);  // other stream

  ModelRandom r(1, 0);
  const Index n = 200000;
  double m = 0.0, v = 0.0;
  for (Index k = 0; k < n; ++k) {
    double x = r.Normal();
    m += x;
    v += x * x;
  }
  m /= double(n);
  v = v / double(n) - m * m;
  BOOST_CHECK_SMALL(m, 0.01);
  BOOST_CHECK_CLOSE(v, 1.0, 1.5);

  for (double mean : {0.5, 30.0, 1700.0}) {
    double s = 0.0, s2 = 0.0;
    const Index draws = 4000;
    for (Index k = 0; k < draws; ++k) {
      double x = double(r.Poisson(mean));
      s += x;
      s2 += x * x;
    }
    double pm = s / double(draws);
    double pv = s2 / double(draws) - pm * pm;
    BOOST_CHECK_SMALL(pm - mean, 5.0 * std::sqrt(mean / double(draws)));
    BOOST_CHECK_CLOSE(pv, mean, 10.0);
  }
  for (Index k = 0; k < 1000; ++k) {
    Index i = r.Integer(7);
    BOOST_CHECK(i >= 0 && i < 7);
  }
}

BOOST_AUTO_TEST_CASE(cells_match_brute_force) {
  // few cells along x (wrapping shells) and many along z
  for (const Eigen::Vector3d& box :
       {Eigen::Vector3d(3.0, 5.0, 11.0), Eigen::Vector3d(2.1, 2.1, 2.1)}) {
    std::vector<Eigen::Vector3d> p = RandomPoints(300, box, 11);
    const double radius = 1.0;
    PeriodicCells cells(box, radius, p);
    for (Index i = 0; i < Index(p.size()); i += 7) {
      std::vector<std::pair<double, Index>> all;
      for (Index j = 0; j < Index(p.size()); ++j) {
        all.emplace_back(MinImage(p[i], p[j], box), j);
      }
      std::sort(all.begin(), all.end());
      for (Index l : {1, 5, 40, 299, 300}) {
        std::vector<Index> nb = cells.Nearest(i, l);
        BOOST_REQUIRE_EQUAL(Index(nb.size()), l);
        BOOST_CHECK_EQUAL(nb[0], i);
        for (Index k = 0; k < l; ++k) {
          BOOST_CHECK_EQUAL(nb[k], all[k].second);
        }
      }
      std::vector<Index> within;
      cells.Within(p[i], 0.8,
                   [&](Index j, const Eigen::Vector3d& conn, double dist) {
                     within.push_back(j);
                     BOOST_CHECK_CLOSE(dist, conn.norm(), 1e-10);
                     BOOST_CHECK_CLOSE(dist, MinImage(p[i], p[j], box), 1e-8);
                   });
      std::sort(within.begin(), within.end());
      std::vector<Index> expected;
      for (const auto& a : all) {
        if (a.first < 0.8) {
          expected.push_back(a.second);
        }
      }
      std::sort(expected.begin(), expected.end());
      BOOST_CHECK(within == expected);
    }
  }
  std::vector<Eigen::Vector3d> p =
      RandomPoints(10, Eigen::Vector3d(4, 4, 4), 1);
  BOOST_CHECK_THROW(PeriodicCells(Eigen::Vector3d(4, 4, 4), 2.0, p),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(point_process_exclusion) {
  StochasticNetwork::Options opt = TestOptions(1500);
  std::vector<Eigen::Vector3d> pos;
  std::vector<double> radii;
  ModelRandom random(1, 0);
  Index it = StochasticNetwork::SamplePoints(opt.points, opt.box, opt.sites,
                                             random, pos, radii);
  BOOST_CHECK(it >= 1);
  BOOST_REQUIRE_EQUAL(Index(pos.size()), opt.sites);
  // no point lies in the sphere of another one
  PeriodicCells cells(opt.box, 1.5 * nm, pos);
  double rmax = *std::max_element(radii.begin(), radii.end());
  BOOST_REQUIRE(rmax < 1.5 * nm);
  Index violations = 0;
  for (Index i = 0; i < opt.sites; ++i) {
    BOOST_CHECK(radii[i] >= opt.points.hard_core);
    cells.Within(pos[i], rmax, [&](Index j, const Eigen::Vector3d&, double d) {
      if (j != i && d < std::max(radii[i], radii[j])) {
        ++violations;
      }
    });
  }
  BOOST_CHECK_EQUAL(violations, 0);

  // too many sites for the box
  StochasticNetwork::PointProcess big = opt.points;
  big.hard_core = 1.0 * nm;
  big.max_iterations = 20;
  ModelRandom r2(1, 0);
  BOOST_CHECK_THROW(
      StochasticNetwork::SamplePoints(big, opt.box, opt.sites, r2, pos, radii),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(moving_average) {
  Eigen::Vector3d box(4.0, 4.0, 4.0);
  std::vector<Eigen::Vector3d> p = RandomPoints(200, box, 2);
  PeriodicCells cells(box, 1.0, p);
  StochasticNetwork::EnergyModel m;
  m.w_a = 0.2;
  m.w_b = 0.5;
  m.l_b = 4;
  m.l_c = 30;
  ModelRandom r(4, 0);
  std::vector<double> A(200), B(200), C(200);
  for (Index i = 0; i < 200; ++i) {
    A[i] = r.Normal();
    B[i] = r.Normal();
    C[i] = r.Normal();
  }
  std::vector<double> z =
      StochasticNetwork::MovingAverage(cells, m, 200, A, B, C);
  for (Index i : {0, 17, 199}) {
    std::vector<Index> nb = cells.Nearest(i, 30);
    double b = 0.0, c = 0.0;
    for (Index k = 0; k < 4; ++k) {
      b += B[nb[k]];
    }
    for (Index k = 0; k < 30; ++k) {
      c += C[nb[k]];
    }
    double expected = std::sqrt(0.2) * A[i] + std::sqrt(0.5 / 4.0) * b +
                      std::sqrt(0.3 / 30.0) * c;
    BOOST_CHECK_CLOSE(z[i], expected, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(network_statistics) {
  StochasticNetwork::Options opt = TestOptions(3000);
  StochasticNetwork net(opt);
  BOOST_REQUIRE_EQUAL(net.NumberOfSites(), 3000);
  // reproducible, and another seed gives another network
  StochasticNetwork same(opt);
  BOOST_CHECK(same.Positions()[123] == net.Positions()[123]);
  BOOST_CHECK_EQUAL(same.SiteEnergies(1)[77], net.SiteEnergies(1)[77]);
  BOOST_REQUIRE_EQUAL(same.Pairs().size(), net.Pairs().size());
  BOOST_CHECK_EQUAL(same.Pairs()[50].J2[0], net.Pairs()[50].J2[0]);
  StochasticNetwork::Options other = opt;
  other.seed = 4;
  BOOST_CHECK(StochasticNetwork(other).Positions()[0] != net.Positions()[0]);

  // site energies: mean, sigma (the moving average has unit variance; the
  // sample sd of a finite network is smaller by about sqrt(1 - cbar))
  for (Index c = 0; c < 2; ++c) {
    const std::vector<double>& e = net.SiteEnergies(c);
    double m = std::accumulate(e.begin(), e.end(), 0.0) / double(e.size());
    double v = 0.0;
    for (double x : e) {
      v += (x - m) * (x - m);
    }
    double sd = std::sqrt(v / double(e.size()));
    const auto& em = opt.carriers[c].energy;
    BOOST_CHECK_SMALL(m - em.mean, 0.3 * em.sigma);
    BOOST_CHECK_CLOSE(sd, em.sigma, 12.0);
  }
  // hole and electron energies are anticorrelated
  {
    const std::vector<double>& a = net.SiteEnergies(0);
    const std::vector<double>& b = net.SiteEnergies(1);
    double ma = std::accumulate(a.begin(), a.end(), 0.0) / 3000.0;
    double mb = std::accumulate(b.begin(), b.end(), 0.0) / 3000.0;
    double sab = 0.0, saa = 0.0, sbb = 0.0;
    for (Index i = 0; i < 3000; ++i) {
      sab += (a[i] - ma) * (b[i] - mb);
      saa += (a[i] - ma) * (a[i] - ma);
      sbb += (b[i] - mb) * (b[i] - mb);
    }
    BOOST_CHECK_CLOSE(sab / std::sqrt(saa * sbb), -0.8, 5.0);
  }

  // pairs: all sites closer than 1 nm are connected, none beyond 2 nm
  Index close = 0;
  PeriodicCells cells(opt.box, 1.0 * nm, net.Positions());
  for (Index i = 0; i < 3000; ++i) {
    cells.Within(
        net.Positions()[i], 1.0 * nm,
        [&](Index j, const Eigen::Vector3d&, double) { close += j > i; });
  }
  Index connected_close = 0;
  double sum_y = 0.0, sum_y2 = 0.0;
  Index n_y = 0;
  for (const auto& p : net.Pairs()) {
    double d = p.R.norm();
    BOOST_CHECK(d < 2.0 * nm);
    BOOST_CHECK(p.i < p.j);
    BOOST_CHECK_CLOSE(
        d, MinImage(net.Positions()[p.i], net.Positions()[p.j], opt.box), 1e-8);
    connected_close += d < 1.0 * nm;
    if (d > 0.9 * nm && d < 1.1 * nm) {
      double y = std::log10(p.J2[0] / (ev * ev));
      sum_y += y;
      sum_y2 += y * y;
      ++n_y;
    }
  }
  BOOST_CHECK_EQUAL(connected_close, close);
  // log10 J^2 near 1 nm: mean -5.3, variance 1.95
  double my = sum_y / double(n_y);
  double vy = sum_y2 / double(n_y) - my * my;
  BOOST_CHECK_SMALL(my - (-5.3), 0.15);
  BOOST_CHECK_CLOSE(vy, 1.95, 15.0);
}

BOOST_AUTO_TEST_CASE(residual_quantiles) {
  StochasticNetwork::CouplingModel m;
  BOOST_CHECK_EQUAL(StochasticNetwork::Residual(m, 0.7), 0.7);
  // a table of normal quantiles reproduces the identity
  const Index K = 200;
  for (Index k = 0; k < K; ++k) {
    m.residual_quantiles.push_back(
        StochasticNetworkFit::NormalQuantile((double(k) + 0.5) / double(K)));
  }
  for (double g : {-4.0, -2.5, -1.0, 0.0, 0.3, 1.7, 3.5}) {
    BOOST_CHECK_SMALL(StochasticNetwork::Residual(m, g) - g, 2e-3);
  }
  // monotone
  double last = -1e9;
  for (double g = -5.0; g < 5.0; g += 0.01) {
    double x = StochasticNetwork::Residual(m, g);
    BOOST_CHECK(x >= last);
    last = x;
  }
  m.residual_quantiles.clear();
  m.mean_offset = -1.0;
  m.mean_slope = -4.0;
  m.var_offset = 1.0;
  m.var_slope = 3.0;
  // (nm2bohr and bohr2nm agree to about 1e-7)
  BOOST_CHECK_CLOSE(StochasticNetwork::Log10J2(m, 1.0 * nm, 0.5),
                    -5.0 + 2.0 * 0.5, 1e-4);
}

BOOST_AUTO_TEST_CASE(topology_and_parameter_file) {
  StochasticNetwork::Options opt = TestOptions(400);
  StochasticNetwork net(opt);
  Topology top;
  net.BuildTopology(top);
  BOOST_REQUIRE_EQUAL(Index(top.Segments().size()), 400);
  BOOST_CHECK_EQUAL(top.NBList().size(), net.Pairs().size());
  const Segment& s = top.Segments()[10];
  BOOST_CHECK_EQUAL(s.getSiteEnergy(QMStateType(QMStateType::Hole)),
                    net.SiteEnergies(0)[10]);
  BOOST_CHECK_EQUAL(s.getSiteEnergy(QMStateType(QMStateType::Electron)),
                    net.SiteEnergies(1)[10]);
  BOOST_CHECK_CLOSE(s.getU_nX_nN(QMStateType(QMStateType::Hole)), 0.12 * ev,
                    1e-10);
  const QMPair* p = top.NBList()[5];
  BOOST_CHECK_EQUAL(p->getJeff2(QMStateType(QMStateType::Electron)),
                    net.Pairs()[5].J2[1]);

  opt.carriers[0].coupling.residual_quantiles = {-2.0, -0.5, 0.1, 0.4, 1.5};
  std::string xml = StochasticNetwork::ToXML(opt);
  const std::string file = "stochasticnetwork_test.xml";
  {
    std::ofstream out(file);
    out << xml;
  }
  votca::tools::Property prop;
  prop.LoadFromXML(file);
  std::remove(file.c_str());
  StochasticNetwork::Options back =
      StochasticNetwork::FromProperty(prop.get("stochasticnetwork"));
  BOOST_CHECK_EQUAL(back.sites, opt.sites);
  BOOST_CHECK_CLOSE(back.box.x(), opt.box.x(), 1e-4);
  BOOST_CHECK_CLOSE(back.points.radius_variance, opt.points.radius_variance,
                    1e-4);
  BOOST_REQUIRE_EQUAL(back.connection_r.size(), opt.connection_r.size());
  BOOST_CHECK_CLOSE(back.connection_r[2], opt.connection_r[2], 1e-4);
  BOOST_REQUIRE_EQUAL(back.carriers.size(), 2);
  BOOST_CHECK(back.carriers[1].state == opt.carriers[1].state);
  BOOST_CHECK_CLOSE(back.carriers[1].energy.sigma, opt.carriers[1].energy.sigma,
                    1e-4);
  BOOST_CHECK_EQUAL(back.carriers[1].energy.l_c, 40);
  BOOST_CHECK_EQUAL(back.carriers[0].coupling.residual_quantiles.size(), 5);
  BOOST_CHECK_CLOSE(back.carriers[0].U_xN_xX, opt.carriers[0].U_xN_xX, 1e-4);
  BOOST_CHECK_CLOSE(back.energy_correlation, -0.8, 1e-4);
  // the network from the file is the same (up to the unit conversions)
  StochasticNetwork again(back);
  BOOST_REQUIRE_EQUAL(again.Pairs().size(), net.Pairs().size());
  BOOST_CHECK_SMALL((again.Positions()[3] - net.Positions()[3]).norm(), 1e-4);
}

BOOST_AUTO_TEST_CASE(fit_weights_exact) {
  StochasticNetworkFit::Curve kappa, ob, oc;
  for (Index k = 0; k < 30; ++k) {
    double r = 0.1 * double(k + 1);
    kappa.r.push_back(r);
    ob.r.push_back(r);
    oc.r.push_back(r);
    ob.value.push_back(std::exp(-r));
    oc.value.push_back(1.0 / (1.0 + r * r));
    kappa.value.push_back(0.3 * ob.value.back() + 0.5 * oc.value.back());
    kappa.count.push_back(double(k + 1));
  }
  double wb = 0.0, wc = 0.0;
  double res = StochasticNetworkFit::FitWeights(kappa, ob, oc, wb, wc);
  // mean square of kappa about 0.1: 1e-14 is far below any real misfit
  BOOST_CHECK_SMALL(res, 1e-14);
  BOOST_CHECK_CLOSE(wb, 0.3, 1e-6);
  BOOST_CHECK_CLOSE(wc, 0.5, 1e-6);
  // constrained: w_b + w_c <= 1
  for (Index k = 0; k < 30; ++k) {
    kappa.value[k] = 0.9 * ob.value[k] + 0.6 * oc.value[k];
  }
  StochasticNetworkFit::FitWeights(kappa, ob, oc, wb, wc);
  BOOST_CHECK(wb >= 0.0 && wc >= 0.0);
  BOOST_CHECK_CLOSE(wb + wc, 1.0, 1e-8);
  // negative correlation: zero weights
  for (Index k = 0; k < 30; ++k) {
    kappa.value[k] = -0.2 * ob.value[k];
  }
  StochasticNetworkFit::FitWeights(kappa, ob, oc, wb, wc);
  BOOST_CHECK_EQUAL(wb, 0.0);
  BOOST_CHECK_EQUAL(wc, 0.0);
}

BOOST_AUTO_TEST_CASE(overlaps) {
  Eigen::Vector3d box(5.0, 5.0, 5.0);
  std::vector<Eigen::Vector3d> p = RandomPoints(250, box, 9);
  std::vector<Index> ls = {1, 5, 50};
  auto o = StochasticNetworkFit::Overlaps(p, box, ls, 0.25, 2.4);
  BOOST_REQUIRE_EQUAL(o.size(), 3);
  // brute force for l = 5 in one bin
  PeriodicCells cells(box, 2.4, p);
  double sum = 0.0, count = 0.0;
  for (Index i = 0; i < 250; ++i) {
    std::vector<Index> a = cells.Nearest(i, 5);
    std::sort(a.begin(), a.end());
    for (Index j = i + 1; j < 250; ++j) {
      double d = MinImage(p[i], p[j], box);
      if (d >= 0.75 && d < 1.0) {
        std::vector<Index> b = cells.Nearest(j, 5);
        std::sort(b.begin(), b.end());
        std::vector<Index> both;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                              std::back_inserter(both));
        sum += double(both.size()) / 5.0;
        count += 1.0;
      }
    }
  }
  BOOST_CHECK_EQUAL(o[1].count[3], count);
  BOOST_CHECK_CLOSE(o[1].value[3], sum / count, 1e-10);
  for (std::size_t b = 0; b < o[0].value.size(); ++b) {
    BOOST_CHECK_EQUAL(o[0].value[b], 0.0);  // only the site itself
    BOOST_CHECK(o[2].value[b] >= 0.0 && o[2].value[b] <= 1.0);
  }
  // larger sets overlap more
  BOOST_CHECK(o[2].value[3] > o[1].value[3]);
}

BOOST_AUTO_TEST_CASE(fit_recovers_parameters) {
  // generate a network with known parameters and fit it again
  StochasticNetwork::Options opt = TestOptions(2500);
  opt.carriers[0].coupling.var_slope = 1.0;
  StochasticNetwork net(opt);
  std::vector<QMStateType> carriers = {QMStateType(QMStateType::Hole),
                                       QMStateType(QMStateType::Electron)};
  StochasticNetworkFit::Data data =
      StochasticNetworkFit::FromNetwork(net, carriers, opt.box);
  StochasticNetworkFit::Settings s;
  s.carriers = carriers;
  s.fit_points = false;
  s.points = opt.points;
  s.l_b = {2, 3, 4, 5, 6, 7, 8, 9, 10};
  for (Index l = 20; l <= 120; l += 10) {
    s.l_c.push_back(l);
  }
  s.correlation_cutoff = 3.0 * nm;
  StochasticNetworkFit fit(data, s);
  StochasticNetworkFit::Result res = fit.Fit();
  const auto& o = res.options;
  // connection probability: 1 below 1 nm, 0.4 at 1.5 nm
  for (std::size_t k = 0; k < o.connection_r.size(); ++k) {
    double r = o.connection_r[k] / nm;
    if (r > 0.6 && r < 0.95) {
      BOOST_CHECK_CLOSE(o.connection_p[k], 1.0, 1e-8);
    }
    if (std::abs(r - 1.525) < 1e-6) {
      BOOST_CHECK_SMALL(o.connection_p[k] - 0.37, 0.05);
    }
  }
  BOOST_CHECK(o.connection_r.back() <= 2.05 * nm);
  // couplings
  const auto& c = o.carriers[0].coupling;
  BOOST_CHECK_SMALL(c.mean_slope - (-5.0), 0.15);
  BOOST_CHECK_SMALL(c.mean_offset - (-0.3), 0.2);
  BOOST_CHECK_SMALL(c.var_slope - 1.0, 0.3);
  BOOST_CHECK_SMALL(o.coupling_correlation - 0.5, 0.05);
  // energies
  BOOST_CHECK_SMALL(o.energy_correlation - (-0.8), 0.05);
  const auto& e = o.carriers[0].energy;
  BOOST_CHECK_CLOSE(e.sigma, 0.09 * ev, 8.0);
  BOOST_CHECK_SMALL(e.w_a - 0.3, 0.08);
  BOOST_CHECK(e.l_b >= 3 && e.l_b <= 10);
  BOOST_CHECK(e.l_c >= 40 && e.l_c <= 90);
  BOOST_CHECK(res.energies[0].rms < 0.03);
  // FromNetwork carries no reorganization energies
  BOOST_CHECK_SMALL(o.carriers[0].U_nX_nN + o.carriers[0].U_xN_xX, 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
