/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE rateengine_test

// Standard includes
#include <cmath>

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/rate_engine.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(rateengine_test)
BOOST_AUTO_TEST_CASE(sign) {

  Segment seg1("null", 0);

  Segment seg2("one", 0);
  Eigen::Vector3d deltaR = {5, 0, 0};
  QMPair pair(0, &seg1, &seg2, deltaR);
  QMStateType e = QMStateType::Electron;
  QMStateType h = QMStateType::Hole;
  seg1.setU_nX_nN(0.002, e);
  seg1.setU_xN_xX(0.001, e);
  seg1.setU_xX_nN(0.001, e);
  seg2.setU_nX_nN(0.002, e);
  seg2.setU_xN_xX(0.001, e);
  seg2.setU_xX_nN(0.001, e);
  seg1.setU_nX_nN(0.002, h);
  seg1.setU_xN_xX(0.001, h);
  seg1.setU_xX_nN(0.001, h);
  seg2.setU_nX_nN(0.002, h);
  seg2.setU_xN_xX(0.001, h);
  seg2.setU_xX_nN(0.001, h);
  pair.setJeff2(1.00e-06, h);
  pair.setJeff2(1.00e-06, e);

  Eigen::Vector3d field = {1.0, 0.0, 0.0};
  field *= 9.72345198649679e-05;
  double temperature = 0.000950043476927;  // 300K
  Rate_Engine engine(temperature, field);
  Rate_Engine::PairRates pr_e = engine.Rate(pair, e);
  Rate_Engine::PairRates pr_h = engine.Rate(pair, h);
  BOOST_CHECK_EQUAL(pr_e.rate12 < pr_e.rate21, true);
  BOOST_CHECK_EQUAL(pr_h.rate12 > pr_h.rate21, true);
}

BOOST_AUTO_TEST_CASE(markus) {

  Segment seg1("null", 0);

  Segment seg2("one", 0);
  Eigen::Vector3d deltaR = {5, 0, 0};
  QMPair pair(0, &seg1, &seg2, deltaR);
  QMStateType e = QMStateType::Electron;
  seg1.setEMpoles(e, -0.003674932248085);
  seg2.setEMpoles(e, -0.005512398372128);
  seg1.setU_nX_nN(0.003674932248085, e);
  seg1.setU_xN_xX(0.004042425472894, e);
  seg1.setU_xX_nN(0.036749322480855, e);
  seg2.setU_nX_nN(0.006614878046554, e);
  seg2.setU_xN_xX(0.004777411922511, e);
  seg2.setU_xX_nN(0.062473848217453, e);
  double J = 1.00e-03;
  pair.setJeff2(J * J, e);

  Eigen::Vector3d field = {1.0, 0.0, 0.0};
  field *= 9.72345198649679e-05;
  double temperature = 0.000950043476927;  // 300K
  Rate_Engine engine(temperature, field);
  Rate_Engine::PairRates pr = engine.Rate(pair, e);
  BOOST_CHECK_CLOSE(pr.rate12, 0.069766184455863, 1e-5);
  BOOST_CHECK_CLOSE(pr.rate21, 221259502589.522, 1e-5);
}

namespace {
// Pair of two segments with identical inner reorganization energies, so that
// reorg12 == reorg21 and the Marcus rates obey detailed balance.
struct SymmetricPair {
  Segment seg1{"one", 0};
  Segment seg2{"two", 1};
  QMPair pair;
  explicit SymmetricPair(const Eigen::Vector3d& dR)
      : pair(0, &seg1, &seg2, dR) {
    for (QMStateType s :
         {QMStateType(QMStateType::Electron), QMStateType(QMStateType::Hole)}) {
      seg1.setU_nX_nN(0.003, s);
      seg1.setU_xN_xX(0.002, s);
      seg2.setU_nX_nN(0.003, s);
      seg2.setU_xN_xX(0.002, s);
      seg1.setEMpoles(s, -0.0040);
      seg2.setEMpoles(s, -0.0047);
      pair.setJeff2(1e-6, s);
    }
  }
};
}  // namespace

BOOST_AUTO_TEST_CASE(outer_sphere_adds_to_both_directions) {
  // lambda_out larger than the inner reorganization energy: subtracting it in
  // one direction would give a negative reorganization energy and NaN rates.
  SymmetricPair p(Eigen::Vector3d(5, 2, 0));
  QMStateType e = QMStateType::Electron;
  p.pair.setLambdaO(0.02, e);
  double temperature = 0.000950043476927;  // 300K
  Rate_Engine engine(temperature, Eigen::Vector3d::Zero());
  Rate_Engine::PairRates pr = engine.Rate(p.pair, e);
  BOOST_CHECK(std::isfinite(pr.rate12) && pr.rate12 > 0.0);
  BOOST_CHECK(std::isfinite(pr.rate21) && pr.rate21 > 0.0);

  // both directions use lambda = inner + outer
  double lambda = 0.003 + 0.002 + 0.02;
  double dG = p.pair.getdE12(e);
  double hbar = votca::tools::conv::hbar * votca::tools::conv::ev2hrt;
  auto marcus = [&](double g) {
    return 2 * votca::tools::conv::Pi / hbar * 1e-6 /
           std::sqrt(4 * votca::tools::conv::Pi * lambda * temperature) *
           std::exp(-(g - lambda) * (g - lambda) / (4 * lambda * temperature));
  };
  BOOST_CHECK_CLOSE(pr.rate12, marcus(dG), 1e-10);
  BOOST_CHECK_CLOSE(pr.rate21, marcus(-dG), 1e-10);
}

BOOST_AUTO_TEST_CASE(detailed_balance) {
  // With reorg12 == reorg21: rate12 / rate21 = exp(dG / kT), where
  // dG = E1 - E2 + q F.R12 is the energy released by the hop 1 -> 2.
  SymmetricPair p(Eigen::Vector3d(5, 2, -1));
  Eigen::Vector3d field(0.7, -0.2, 0.4);
  field *= 9.72345198649679e-05;
  double temperature = 0.000950043476927;  // 300K
  Rate_Engine engine(temperature, field);
  for (double lambda_out : {0.0, 0.004, 0.02}) {
    for (QMStateType s :
         {QMStateType(QMStateType::Electron), QMStateType(QMStateType::Hole)}) {
      p.pair.setLambdaO(lambda_out, s);
      double charge = (s == QMStateType::Electron) ? -1.0 : 1.0;
      double dG = p.pair.getdE12(s) + charge * p.pair.R().dot(field);
      Rate_Engine::PairRates pr = engine.Rate(p.pair, s);
      BOOST_CHECK_CLOSE(std::log(pr.rate12 / pr.rate21), dG / temperature,
                        1e-8);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
