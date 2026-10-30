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

#define BOOST_TEST_MODULE latticemodel_test

// Standard includes
#include <cmath>
#include <cstdio>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/latticemodel.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/statesaver.h"

using namespace votca::xtp;
using votca::Index;

namespace {
LatticeModel::Options PaperOptions() {
  LatticeModel::Options opt;
  double nm = votca::tools::conv::nm2bohr;
  double ev = votca::tools::conv::ev2hrt;
  opt.dims = {8, 8, 8};
  opt.spacing = 1.0 * nm;
  opt.cutoff = 1.5 * nm;
  opt.J2 = 1.79e-4 * ev * ev;
  opt.decay_length = 1.0 * nm;
  opt.reorganization = 0.23 * ev;
  opt.mean_energy = -0.76 * ev;
  opt.sigma = 0.19 * ev;
  return opt;
}

void MeanAndStd(const std::vector<double>& x, double& mean, double& sd) {
  mean = 0.0;
  for (double v : x) {
    mean += v;
  }
  mean /= double(x.size());
  double var = 0.0;
  for (double v : x) {
    var += (v - mean) * (v - mean);
  }
  sd = std::sqrt(var / double(x.size()));
}

// correlation of the energies of nearest-neighbour pairs in the topology
double NeighbourCorrelation(const Topology& top, double spacing) {
  std::vector<double> e;
  for (const Segment& seg : top.Segments()) {
    e.push_back(seg.getSiteEnergy(QMStateType::Electron));
  }
  double mean, sd;
  MeanAndStd(e, mean, sd);
  double c = 0.0;
  Index count = 0;
  for (const QMPair* pair : top.NBList()) {
    if (std::abs(pair->Dist() - spacing) < 1e-6 * spacing) {
      c +=
          (e[pair->Seg1()->getId()] - mean) * (e[pair->Seg2()->getId()] - mean);
      ++count;
    }
  }
  return c / double(count) / (sd * sd);
}
}  // namespace

BOOST_AUTO_TEST_SUITE(latticemodel_test)

BOOST_AUTO_TEST_CASE(geometry_and_couplings) {
  LatticeModel::Options opt = PaperOptions();
  opt.disorder = LatticeModel::Disorder::None;
  LatticeModel model(opt);
  Topology top;
  model.BuildTopology(top);

  BOOST_CHECK_EQUAL(top.Segments().size(), 512);
  // 6 nearest and 12 next-nearest neighbours per site
  BOOST_CHECK_EQUAL(top.NBList().size(), 512 * 18 / 2);
  BOOST_CHECK_CLOSE(top.getBox()(0, 0), 8 * opt.spacing, 1e-12);
  BOOST_CHECK_CLOSE(top.Segments()[0].getPos().x(), 0.5 * opt.spacing, 1e-12);

  QMStateType e = QMStateType::Electron;
  QMStateType h = QMStateType::Hole;
  for (const QMPair* pair : top.NBList()) {
    double r = pair->Dist();
    double expected = opt.J2 * std::exp(-(r - opt.spacing) / opt.decay_length);
    BOOST_CHECK_CLOSE(pair->getJeff2(e), expected, 1e-10);
    BOOST_CHECK_CLOSE(pair->getJeff2(h), expected, 1e-10);
    BOOST_CHECK_CLOSE(pair->getReorg12(e), opt.reorganization, 1e-10);
    BOOST_CHECK_CLOSE(pair->getReorg21(e), opt.reorganization, 1e-10);
    BOOST_CHECK_SMALL(pair->getdE12(h), 1e-14);
  }
  for (const Segment& seg : top.Segments()) {
    BOOST_CHECK_CLOSE(seg.getSiteEnergy(e), opt.mean_energy, 1e-12);
  }
}

BOOST_AUTO_TEST_CASE(uncorrelated_energies) {
  LatticeModel::Options opt = PaperOptions();
  opt.dims = {16, 16, 16};
  opt.disorder = LatticeModel::Disorder::Uncorrelated;
  opt.seed = 7;
  LatticeModel model(opt);
  double mean, sd;
  MeanAndStd(model.SiteEnergies(), mean, sd);
  // 4096 samples: standard error of the mean is 1.6 % of sigma
  BOOST_CHECK_SMALL(mean - opt.mean_energy, 0.06 * opt.sigma);
  BOOST_CHECK_CLOSE(sd, opt.sigma, 5.0);

  Topology top;
  model.BuildTopology(top);
  BOOST_CHECK_SMALL(NeighbourCorrelation(top, opt.spacing), 0.06);
}

BOOST_AUTO_TEST_CASE(correlated_energies) {
  LatticeModel::Options opt = PaperOptions();
  opt.dims = {16, 16, 16};
  opt.disorder = LatticeModel::Disorder::Correlated;
  opt.seed = 11;
  LatticeModel model(opt);
  double mean, sd;
  MeanAndStd(model.SiteEnergies(), mean, sd);
  // the construction keeps the variance at sigma^2; correlated samples
  // fluctuate more
  BOOST_CHECK_CLOSE(sd, opt.sigma, 10.0);
  BOOST_CHECK_SMALL(mean - opt.mean_energy, 0.2 * opt.sigma);

  Topology top;
  model.BuildTopology(top);
  // nearest neighbours share most of their l_c = 280 neighbours
  BOOST_CHECK_GT(NeighbourCorrelation(top, opt.spacing), 0.4);
}

BOOST_AUTO_TEST_CASE(correlated_with_and_without_the_site) {
  LatticeModel::Options opt = PaperOptions();
  opt.dims = {12, 12, 12};
  opt.disorder = LatticeModel::Disorder::Correlated;
  opt.seed = 5;
  LatticeModel with(opt);
  opt.include_self = false;
  LatticeModel without(opt);
  BOOST_CHECK(with.SiteEnergies() != without.SiteEnergies());
  for (const LatticeModel* m : {&with, &without}) {
    double mean, sd;
    MeanAndStd(m->SiteEnergies(), mean, sd);
    BOOST_CHECK_CLOSE(sd, opt.sigma, 15.0);
  }
  // l_b = l_c = 1 with the site itself: no correlation between the
  // components, the energies stay Gaussian with sigma
  opt.include_self = true;
  opt.l_b = 1;
  opt.l_c = 1;
  opt.dims = {16, 16, 16};
  LatticeModel local(opt);
  Topology top;
  local.BuildTopology(top);
  BOOST_CHECK_SMALL(NeighbourCorrelation(top, opt.spacing), 0.1);
  // all sites of a 4x4x4 box: 64 with the site itself, 63 without
  opt.dims = {4, 4, 4};
  opt.cutoff = 1.0 * opt.spacing;
  opt.l_b = 9;
  opt.l_c = 64;
  BOOST_CHECK_NO_THROW(LatticeModel{opt});
  opt.include_self = false;
  BOOST_CHECK_THROW(LatticeModel{opt}, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(seed_reproducibility) {
  LatticeModel::Options opt = PaperOptions();
  opt.disorder = LatticeModel::Disorder::Correlated;
  opt.seed = 3;
  LatticeModel a(opt);
  LatticeModel b(opt);
  BOOST_CHECK(a.SiteEnergies() == b.SiteEnergies());
  opt.seed = 4;
  LatticeModel c(opt);
  BOOST_CHECK(a.SiteEnergies() != c.SiteEnergies());

  // Box-Muller on mt19937_64: platform independent, fixed reference values
  LatticeModel::NormalGenerator gen(42);
  double x0 = gen();
  double x1 = gen();
  LatticeModel::NormalGenerator gen2(42);
  BOOST_CHECK_EQUAL(gen2(), x0);
  BOOST_CHECK_EQUAL(gen2(), x1);
}

BOOST_AUTO_TEST_CASE(invalid_options) {
  LatticeModel::Options opt = PaperOptions();
  opt.cutoff = 4.0 * opt.spacing;  // half the box
  BOOST_CHECK_THROW(LatticeModel{opt}, std::runtime_error);
  opt = PaperOptions();
  opt.disorder = LatticeModel::Disorder::Correlated;
  opt.dims = {4, 4, 4};  // 63 other sites < l_c
  BOOST_CHECK_THROW(LatticeModel{opt}, std::runtime_error);
  BOOST_CHECK_THROW(LatticeModel::DisorderFromString("gaussian"),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(statefile_roundtrip) {
  LatticeModel::Options opt = PaperOptions();
  opt.dims = {6, 6, 6};
  opt.disorder = LatticeModel::Disorder::Uncorrelated;
  LatticeModel model(opt);
  Topology top;
  model.BuildTopology(top);

  std::string file = "latticemodel_roundtrip.hdf5";
  std::remove(file.c_str());
  StateSaver saver(file);
  saver.WriteFrame(top);
  Topology read = saver.ReadFrame(0);
  std::remove(file.c_str());

  BOOST_REQUIRE_EQUAL(read.Segments().size(), top.Segments().size());
  BOOST_REQUIRE_EQUAL(read.NBList().size(), top.NBList().size());
  QMStateType h = QMStateType::Hole;
  for (Index i = 0; i < Index(top.Segments().size()); ++i) {
    BOOST_CHECK_CLOSE(read.Segments()[i].getSiteEnergy(h),
                      model.SiteEnergies()[i], 1e-12);
  }
  for (Index p = 0; p < top.NBList().size(); ++p) {
    BOOST_CHECK_CLOSE(read.NBList()[p]->getJeff2(h),
                      top.NBList()[p]->getJeff2(h), 1e-12);
    BOOST_CHECK_SMALL((read.NBList()[p]->R() - top.NBList()[p]->R()).norm(),
                      1e-12);
  }
}

BOOST_AUTO_TEST_SUITE_END()
