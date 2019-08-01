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
BOOST_AUTO_TEST_CASE(markus) {

  Segment seg1("null", 0);

  Segment seg2("one", 0);
  Eigen::Vector3d deltaR = {2, 1.5, 3};
  QMPair pair(0, &seg1, &seg2, deltaR);
  QMStateType e = QMStateType::Electron;
  seg1.setEMpoles(e, 0.3);
  seg2.setEMpoles(e, 0.4);
  seg1.setU_nX_nN(0.05, e);
  seg1.setU_xN_xX(0.02, e);
  seg1.setU_xX_nN(0.8, e);
  seg2.setU_nX_nN(0.04, e);
  seg2.setU_xN_xX(0.02, e);
  seg2.setU_xX_nN(0.7, e);
  pair.setJeff2(0.1, e);

  Eigen::Vector3d field = {1.0, 1.0, 1.0};
  double temperature =
      300 * votca::tools::conv::kB * votca::tools::conv::ev2hrt;
  Rate_Engine engine(temperature, field);
  Rate_Engine::PairRates pr = engine.Rate(pair, e);
  std::cout << pr.rate12 << std::endl;
  std::cout << pr.rate21 << std::endl;
  BOOST_CHECK_CLOSE(pr.rate12, 0.0, 1e-6);
  BOOST_CHECK_CLOSE(pr.rate21, 0.0, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
