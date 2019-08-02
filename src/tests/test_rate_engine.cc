/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <votca/xtp/rate_engine.h>

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

BOOST_AUTO_TEST_SUITE_END()
