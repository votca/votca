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

#define BOOST_TEST_MODULE qmpair_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/qmpair.h>

using namespace votca::tools;

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(qmpair_test)

BOOST_AUTO_TEST_CASE(getters_test) {

  Segment seg("one", 1);
  Atom atm1(1, "C", Eigen::Vector3d::Ones());
  seg.push_back(atm1);

  Segment seg2("two", 2);
  Atom atm2(2, "C", -Eigen::Vector3d::Ones());
  seg2.push_back(atm2);

  QMPair pair(0, &seg, &seg2, 0.5 * Eigen::Vector3d::Ones());

  pair.setType(QMPair::PairType::Hopping);

  BOOST_CHECK_EQUAL(pair.getType(), QMPair::Hopping);
}

BOOST_AUTO_TEST_SUITE_END()
