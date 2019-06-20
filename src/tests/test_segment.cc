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

#define BOOST_TEST_MODULE segment_test
#include <boost/test/unit_test.hpp>
#include <vector>
#include <votca/xtp/segment.h>

using namespace votca::tools;
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(segment_test)

BOOST_AUTO_TEST_CASE(getElementtest) {
  Segment seg("one", 1);
  Atom atm1(1, "C", Eigen::Vector3d::Zero());
  Atom atm2(2, "H", Eigen::Vector3d::UnitX());
  Atom atm3(3, "N", Eigen::Vector3d::UnitY());
  Atom atm4(4, "Au", Eigen::Vector3d::UnitZ());
  Atom atm5(5, "C", 2 * Eigen::Vector3d::UnitZ());
  seg.FindUniqueElements();
}

BOOST_AUTO_TEST_SUITE_END()
